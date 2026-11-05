# HybridPIR

Implementation of the Minimal and Compressed configurations of HybridPIR
(p = 2^4, n = 2^20 records of 256 B, N = 2048, q = 2^64).

## Requirements

- Linux with `g++` (OpenMP) and `taskset`
- x86-64 CPU with AVX2 and FMA, at least 4 cores
- about 7 GiB of free memory

## Build

```
make
```

## Run

```
./run.sh minimal 1       # Minimal,    1 thread
./run.sh minimal 4       # Minimal,    4 threads
./run.sh compressed 1    # Compressed, 1 thread
./run.sh compressed 4    # Compressed, 4 threads
```

Minimal applies record packing, FFT preprocessing and modulus switching of the
response. Compressed additionally sends the query mask as a PRG seed and
ring-switches the response to a ring of dimension 1024.

Each run generates a random database and fresh keys, queries a random index,
and ends with a `SUMMARY` block:

| Field | Meaning |
|---|---|
| `Query(total)` | query size sent by the client |
| `EvalKey size (coef)` | evaluation keys sent by the client |
| `Response after  RS` | response size |
| `preprocessing_ms` | database encoding and FFT |
| `query_unpack_ms` | query unpacking |
| `main_comp_ms` | first dimension and CMUX tree |
| `extract_time_ms` | packed-record extraction and response compression |
| `total_ms` | online latency |
| `record match?` | whether the decrypted record equals the requested one |

`spqlios/` is taken from the TFHE library.
