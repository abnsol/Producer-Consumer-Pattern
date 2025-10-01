# Producer-Consumer Refactor

## Overview

In `analysis_tasks.py`, the fine-mapping batch worker (`finemap_region_batch_worker`) was refactored to use a producer-consumer threading pattern. Previously, the worker both computed fine-mapping results and saved them to the database in the same thread. Now, the computation (producer) and database saving (consumer) are separated into two threads, communicating via a shared queue.

- **Decouple CPU and I/O:** Fine-mapping is CPU-intensive, while saving to the database is I/O-bound. By separating these tasks, CPU computation is not blocked by slow I/O operations.
- **Maintainability:** The producer-consumer pattern is a well-established design for handling such workflows, making the codebase easier to understand and extend.
- **No Nested Processes:** Threads are used instead of processes to avoid the complexity and overhead of nested multiprocessing.

## How It Works

- **Producer Thread:** Runs the fine-mapping (`finemap_region`) for each region and puts the result (a DataFrame of credible SNPs) onto a queue.
- **Consumer Thread:** Pulls results from the queue and saves them to the database using `db.save_lead_variant_credible_sets`.
- **Batch Results:** Results are still accumulated in memory (`batch_results`) and returned at the end, as downstream tasks require the full batch output.

## Considerations

- **Memory Usage:** The memory footprint is similar to the previous approach, as all credible SNPs for the batch are still stored in memory. The producer-consumer pattern does not reduce memory usage but improves efficiency by decoupling computation from I/O.
- **No Memory Bloat Expected:** The number of credible SNPs per region is typically small, so memory usage should remain manageable unless processing extremely large batches.
- **No Nested Processes:** Using threads avoids the complexity of managing multiple process pools.

## Resources

- [Producer-Consumer Pattern Example Blog](https://lionlai1989.github.io/python/multiprocessing/ProcessQueue-and-Producer-Consumer-Pattern)
