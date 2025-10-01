import asyncio
import os
import time
from flask import json
from loguru import logger
from prefect import flow
from status_tracker import TaskState
import multiprocessing as mp

from analysis_tasks import (
    munge_sumstats_preprocessing, filter_significant_variants, run_cojo_per_chromosome, create_region_batches, finemap_region_batch_worker,
    save_sumstats_for_workers, cleanup_sumstats_file
)

from project_tasks import (
    save_analysis_state_task, create_analysis_result_task, 
    get_project_analysis_path_task
)

import pandas as pd
from datetime import datetime, timezone
from prefect.task_runners import ThreadPoolTaskRunner

from config import Config

@flow(log_prints=True)
def analysis_pipeline_flow(projects_handler, analysis_handler, mongodb_uri, db_name, user_id, project_id, gwas_file_path, ref_genome="GRCh37", 
                           population="EUR", batch_size=5, max_workers=3,
                           maf_threshold=0.01, seed=42, window=2000, L=-1, 
                           coverage=0.95, min_abs_corr=0.5):
    """
    Complete analysis pipeline flow using Prefect for orchestration
    but multiprocessing for fine-mapping batches (R safety)
    """
    
    logger.info(f"[PIPELINE] Starting Prefect analysis pipeline with multiprocessing fine-mapping")
    logger.info(f"[PIPELINE] Project: {project_id}, User: {user_id}")
    logger.info(f"[PIPELINE] File: {gwas_file_path}")
    logger.info(f"[PIPELINE] Batch size: {batch_size} regions per worker process")
    logger.info(f"[PIPELINE] Max workers: {max_workers}")
    logger.info(f"[PIPELINE] Parameters: maf={maf_threshold}, seed={seed}, window={window}kb, L={L}, coverage={coverage}, min_abs_corr={min_abs_corr}")
    
    try:
        # Get project-specific output directory (using Prefect task)
        output_dir = get_project_analysis_path_task.submit(projects_handler, user_id, project_id).result()
        logger.info(f"[PIPELINE] Using output directory: {output_dir}")
        
        # Save initial analysis state
        initial_state = {
            "status": "Running",
            "stage": "Preprocessing",
            "progress": 10,
            "message": "Starting MungeSumstats preprocessing",
            "started_at": datetime.now(timezone.utc).isoformat(),
        }
        save_analysis_state_task.submit(projects_handler, user_id, project_id, initial_state).result()
        
        logger.info(f"[PIPELINE] Stage 1: MungeSumstats preprocessing")
        # munged_file_result = munge_sumstats_preprocessing.submit(gwas_file_path, output_dir, ref_genome=ref_genome, n_threads=14).result()
        from scripts.mock_munge import mock_format_sumstats
        # inside your munge_sumstats_preprocessing task or before calling the flow:
        munged_df, munged_file_result = mock_format_sumstats(
            input_path="data/gwas_sumstats/mock_gwas.tsv",
            output_dir="/workspaces/Producer-Consumer-Pattern/local_runs/user_demo_user/project_demo_proj",
            output_name="munged_sumstats.tsv",
        )
        
        # Extract the actual file path from the result
        if isinstance(munged_file_result, tuple):
            munged_df, munged_file = munged_file_result
        else:
            munged_file = munged_file_result
            munged_df = pd.read_csv(munged_file, sep='\t')
        
        # Update analysis state after preprocessing
        preprocessing_state = {
            "status": "Running",
            "stage": "Filtering",
            "progress": 30,
            "message": "Preprocessing completed, filtering significant variants",
            "started_at": initial_state["started_at"]
        }
        save_analysis_state_task.submit(projects_handler, user_id, project_id, preprocessing_state).result()
        
        logger.info(f"[PIPELINE] Stage 2: Loading and filtering variants")
        significant_df_result = filter_significant_variants.submit(munged_df, output_dir).result()
        
        # Extract the actual DataFrame
        if isinstance(significant_df_result, tuple):
            significant_df, _ = significant_df_result
        else:
            significant_df = significant_df_result
        
        # Update analysis state after filtering
        filtering_state = {
            "status": "Running",
            "stage": "Cojo",
            "progress": 50,
            "message": "Filtering completed, running COJO analysis"
        }
        save_analysis_state_task.submit(projects_handler, user_id, project_id, filtering_state).result()
        
        logger.info(f"[PIPELINE] Stage 3: COJO analysis")
       
        config = Config.from_env()
        plink_dir = config.plink_dir
        cojo_result = run_cojo_per_chromosome.submit(significant_df, plink_dir, output_dir, maf_threshold=maf_threshold, population=population).result()
        
        # Extract the actual DataFrame
        if isinstance(cojo_result, tuple):
            cojo_results, _ = cojo_result
        else:
            cojo_results = cojo_result
        
        if cojo_results is None or len(cojo_results) == 0:
            logger.error("[PIPELINE] No COJO results to process")
            # Save failed state
            failed_state = {
                "status": "Failed",
                "stage": "Cojo",
                "progress": 50,
                "message": "COJO analysis failed - no independent signals found",
            }
            save_analysis_state_task.submit(projects_handler, user_id, project_id, failed_state).result()
            return None
        
        # Update analysis state after COJO
        cojo_state = {
            "status": "Running",
            "stage": "Fine_mapping",
            "progress": 70,
            "message": "COJO analysis completed, starting fine-mapping"
        }
        save_analysis_state_task.submit(projects_handler, user_id, project_id, cojo_state).result()
        
        logger.info(f"[PIPELINE] Stage 4: Multiprocessing fine-mapping)")
        logger.info(f"[PIPELINE] Processing {len(cojo_results)} regions with {batch_size} regions per batch")
        region_batches = create_region_batches(cojo_results, batch_size=batch_size)
        logger.info(f"[PIPELINE] Created {len(region_batches)} batches for {max_workers} worker processes")
        
        sumstats_temp_file = save_sumstats_for_workers(significant_df, output_dir)
        
        # Prepare batch data for multiprocessing
        batch_data_list = []
        for i, batch in enumerate(region_batches):
            db_params = {
                'uri': mongodb_uri,
                'db_name': db_name
            }
            batch_data = (batch, f"batch_{i}", sumstats_temp_file, {
                'db_params': db_params,
                'user_id': user_id,
                'project_id': project_id,
                'finemap_params': {
                    'seed': seed,
                    'window': window,
                    'L': L,
                    'coverage': coverage,
                    'min_abs_corr': min_abs_corr,
                    'population': population
                }
            })
            batch_data_list.append(batch_data)
        
        original_method = mp.get_start_method()
        if original_method != 'spawn':
            logger.info(f"[PIPELINE] Switching multiprocessing method from '{original_method}' to 'spawn' to reduce memory usage")
            mp.set_start_method('spawn', force=True)
        
        all_results = []
        successful_batches = 0
        failed_batches = 0
        
        try:
            with mp.Pool(max_workers) as pool:
                try:
                    # Process all batches in parallel
                    batch_results_list = pool.map(finemap_region_batch_worker, batch_data_list)
                    
                    # Collect results
                    for i, batch_results in enumerate(batch_results_list):
                        if batch_results and len(batch_results) > 0:
                            all_results.extend(batch_results)
                            successful_batches += 1
                            logger.info(f"[PIPELINE] Batch {i} completed with {len(batch_results)} regions")
                        else:
                            failed_batches += 1
                            logger.warning(f"[PIPELINE] Batch {i} failed or returned no results")
                            
                except Exception as e:
                    logger.error(f"[PIPELINE] Error in multiprocessing: {str(e)}")
                    raise
                finally:
                    # Clean up temporary sumstats file after all workers are done
                    cleanup_sumstats_file(sumstats_temp_file)
        finally:
            # Restore original multiprocessing method
            if original_method != 'spawn':
                try:
                    mp.set_start_method(original_method, force=True)
                    logger.info(f"[PIPELINE] Restored multiprocessing method to '{original_method}'")
                except:
                    logger.warning(f"[PIPELINE] Could not restore multiprocessing method to '{original_method}'")
        
        # Combine and save results
        if all_results:
            logger.info(f"[PIPELINE] Combining results from {successful_batches} successful batches")
            combined_results = pd.concat(all_results, ignore_index=True)
            
            # Save results using Prefect tasks
            results_file = create_analysis_result_task.submit(analysis_handler, user_id, project_id, combined_results, output_dir).result()
            
            # Summary statistics
            total_variants = len(combined_results)
            high_pip_variants = len(combined_results[combined_results['PIP'] > 0.5])
            total_credible_sets = combined_results.get('credible_set', pd.Series([0])).max()
            
            # Save completed analysis state
            completed_state = {
                "status": "Completed",
                "progress": 100,
                "message": "Analysis completed successfully",
            }
            save_analysis_state_task.submit(projects_handler, user_id, project_id, completed_state).result()
            
            logger.info(f"[PIPELINE] Analysis completed successfully!")
            logger.info(f"[PIPELINE] - Total variants: {total_variants}")
            logger.info(f"[PIPELINE] - High-confidence variants (PIP > 0.5): {high_pip_variants}")
            logger.info(f"[PIPELINE] - Total credible sets: {total_credible_sets}")
            logger.info(f"[PIPELINE] - Successful batches: {successful_batches}/{len(region_batches)}")
            logger.info(f"[PIPELINE] - Results saved: {results_file}")
            
            return {
                "results_file": results_file,
                "total_variants": total_variants,
                "high_pip_variants": high_pip_variants,
                "total_credible_sets": total_credible_sets
            }
        else:
            logger.error("[PIPELINE]  No fine-mapping results generated")
            # Save failed state for fine-mapping
            failed_finemap_state = {
                "status": "Failed",
                "stage": "Fine_mapping",
                "progress": 70,
                "message": "Fine-mapping failed - no results generated",
            }
            save_analysis_state_task.submit(projects_handler, user_id, project_id, failed_finemap_state).result()
            raise RuntimeError("All fine-mapping batches failed")
            
    except Exception as e:
        logger.error(f"[PIPELINE]  Analysis pipeline failed: {str(e)}")
        # Save failed analysis state
        try:
            failed_state = {
                "status": "Failed",
                "stage": "Unknown",
                "progress": 0,
                "message": f"Analysis pipeline failed: {str(e)}",
            }
            save_analysis_state_task.submit(projects_handler, user_id, project_id, failed_state).result()
        except Exception as state_e:
            logger.error(f"[PIPELINE] Failed to save error state: {str(state_e)}")
        raise