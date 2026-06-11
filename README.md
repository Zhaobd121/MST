# MST (Set) + MST-vector (Vector)

This repository provides two C++ implementations of an MST indexing and querying pipeline:

- **MST** (folder: `MST/`) — designed for **set datasets**.  
  In this version, **epsilon is provided as a fraction**: `<epsilon_num> <epsilon_denom>`.

- **MST-vector** (folder: `MST-vector/`) — designed for **vector datasets**.  
  In this version, **epsilon is provided as a floating-point number**: `epsilon`.

Both versions generally follow the same workflow:
1) **Index** the dataset  
2) **Query** using either `exact` or `density`

---

## Prerequisites

- A C++ compiler, e.g., `g++`
- `make`

> If you encounter a stack overflow due to deep recursion, run:
> ```bash
> ulimit -s unlimited
> ```

---

## Project Layout

- `MST/` — set dataset version
- `MST-vector/` — vector dataset version

---

# Part 1 — MST Set Dataset

## Build MST

Build MST in its own directory.

```bash
cd MST
make
```

The executable name depends on your Makefile, often `program`.

If yours is different, replace `./program` in the commands below.

## Input Format

Example format:

```txt
<n1> <n2> <m>
<left_object_1>:<right_1> <right_2> ...
<left_object_2>:<right_1> <right_5> ...
...
```

The first line provides dataset meta information.

Each subsequent line describes a left object and its associated right objects.

## Usage

### 1) Build index

```bash
./program index <miu> <input_file>
```

### 2) Query

```bash
./program query <exact|density> <miu> <epsilon_num> <epsilon_denom> <input_file>
```

- `miu`: algorithm threshold
- `epsilon_num epsilon_denom`: epsilon as a fraction

Example: `1 10` means `epsilon = 0.1`.

---

# Part 2 — MST-vector Vector Dataset

## Build MST-vector

Build MST-vector in its own directory.

```bash
cd MST-vector
make
```

The executable name depends on your Makefile, often `program`.

If yours is different, replace `./program` in the commands below.

## Input Format

Example format:

```txt
<n>
<vector_0>:<value_0_1> <value_0_2> ...
<vector_1>:<value_1_1> <value_1_2> ...
...
```

The first line is the number of vectors.

Each following line contains a vector identifier and its values.

## Usage

### 1) Build index

```bash
./program index <miu> <input_file>
```

### 2) Query

```bash
./program query <exact|density> <miu> <epsilon> <input_file>
```

- `epsilon`: a floating-point value, e.g., `0.1`

---

## Output

The output format depends on the mode, `exact` or `density`, and the implementation.

Please refer to the code and comments for details if you need to parse results programmatically.

---

## Notes

If you want GitHub to display this documentation on the repository homepage, keep this file at the repository root as `README.md`.

Keep subfolder READMEs if you want module-level docs; otherwise, this integrated README is enough.
