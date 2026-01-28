# Clustering Index and Query Tool

This C++ tool builds indexes and performs set clustering (either **exact** or **density-based**) using a fast, scalable algorithm.

---

## 🔧 Requirements

- C++17 compatible compiler
- Linux or macOS environment
- Dataset in the correct format (see below)

---

## ⚙️ Compilation

Use the provided Makefile:

```bash
make clean && make
```

---

## 📁 Dataset Format

The dataset represents a vector dataset. It must follow this format:

```
<n>
<vector_0>:<value0_1> <value0_2> ...
<vector_1>:<value1_1> <value1_2> ...
...
```

- `n`: number of vectors in the dataset
- vector_num **must be numbered from `0` to `n-1`** without gaps

**Example:**

```
13
0: 8 1
1: 11 1
2: 9 5
3: 7 5
4: 1 5
5: 0 1
6: 0 5
7: 2 1
8: 1 0
9: 13 5
10: 4 6
11: 7 0
12: 14 7
```

---

## 🚀 Usage

### 🏗️ Indexing Mode

Builds the binary core index files.

```bash
./program index <miu> <input_file>
```

**Arguments:**
- `<miu>`: size threshold used for index construction
- `<input_file>`: path to the dataset file

**Output:**

```
index/<input_file_stem>_<miu>/
├── newFCI.txt
├── newNCI.txt
└── newDNCI.txt
```

Each text file stores edges as:  
`(u:int32, v:int32, num:int32)`

---

### 🔍 Query Mode

Runs clustering using either exact or density-based methods:

```bash
./program query <exact|density> <miu> <epsilon> <input_file>
```

**Arguments:**
- `<exact|density>`: query type
- `<miu>`: size threshold used in the index
- `<epsilon>`: distance threshold used in the index
- `<input_file>`: path to dataset

**Output:**

Results are saved to:

```
result/<input_file_stem>_<method>_<miu>_<epsilon>
```

---

## 📄 Output Format

Each cluster is printed like this:

```
Cluster rep: 123
  Cores: 123 124 125
  Non-Cores: 200 201
```

---

## ⚠️ Notes

- You **must run index mode first** before querying.
- If your dataset is large, the algorithm may hit the system stack size limit due to recursion when building the KRT.
- To avoid stack overflow, run this **before execution**:

```bash
ulimit -s unlimited
```

---

## 🧪 Example

```bash
./program index 5 dataset/toyset.txt
./program index 4 dataset/toyset.txt
./program query exact 5 4.8 dataset/toyset.txt
./program query exact 3 5 dataset/toyset.txt
```