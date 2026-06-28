# EGFR Exploratory Custom PyTorch GCN Baseline Report

An exploratory custom PyTorch dense GCN baseline was run using padded RDKit molecular graphs. It is retained as negative benchmark evidence because it did not outperform the Morgan Random Forest baseline.

- GNN status: completed
- Backend used: custom_pytorch_dense_gcn
- Torch Python: `python`
- Torch version: 2.11.0+cu128
- CUDA available to Torch: True
- Device used: NVIDIA GeForce RTX 4090
- Graph rows: 10,593
- Max atoms: 64
- Node feature dimension: 15

## Baseline Metrics

| split | MAE | RMSE | R2 | Pearson | Spearman |
| --- | --- | --- | --- | --- | --- |
| random_split | 0.886 | 1.115 | 0.310 | 0.575 | 0.573 |
| scaffold_split | 0.909 | 1.149 | 0.198 | 0.479 | 0.463 |

## Morgan Random Forest References

- Random split RF: MAE 0.516, RMSE 0.712, R2 0.719
- Scaffold split RF: MAE 0.667, RMSE 0.871, R2 0.550

- GNN beat Morgan RF on random RMSE: False
- GNN beat Morgan RF on scaffold RMSE: False

This is not presented as a mature graph-model result.
