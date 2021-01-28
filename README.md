# Removing many hashes from a truckload of queries

## Running

```
cargo run --release -- \
  -k 57 \
  --scaled 100 \
  data/GCA_001593925.1_ASM159392v1_protein.faa.gz.sig \
  <(find data -type f)
```
