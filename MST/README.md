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

The dataset represents a set. It must follow this format:

```
<n1> <n2> <m>
<left_object_1>:<right_1> <right_2> ...
<left_object_2>:<right_1> <right_2> ...
...
```

- `n1`: number of unique objects on the **left**
- `n2`: number of unique objects on the **right**
- `m`: total number of pairs (connections)
- Each of the next `n1` lines describes a left object’s connections
- Left objects **must be numbered from `1` to `n1`** without gaps

**Example:**

```
5 6 8
1:10
2:11 12
3:13 14 15
4:16
5:17
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
├── FCI.bin
├── NCI.bin
└── DNCI.bin
```

Each binary file stores edges as 16 bytes:  
`(u:int32, v:int32, num:int32, denom:int32)`

---

### 🔍 Query Mode

Runs clustering using either exact or density-based methods:

```bash
./program query <exact|density> <miu> <epsilon_num> <epsilon_denom> <input_file>
```

**Arguments:**
- `<exact|density>`: query type
- `<miu>`: size threshold used in the index
- `<epsilon_num>` and `<epsilon_denom>`: threshold as a fraction ε = num / denom
- `<input_file>`: path to dataset

**Output:**

Results are saved to:

```
result/<input_file_stem>_<method>_<miu>_<epsilon_num>_<epsilon_denom>
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
./program index 3 dataset/actor_movie_set.txt
./program index 7 dataset/actor_movie_set.txt
./program query exact 3 2 3 dataset/actor_movie_set.txt
./program query exact 7 1 2 dataset/actor_movie_set.txt
```