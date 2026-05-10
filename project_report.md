# HelixForge — Genome Assembler
## Project Report | Design and Analysis of Algorithms (DAA)

---

> **Course:** Design and Analysis of Algorithms  
> **Project Title:** HelixForge — A De Bruijn Graph-Based Genome Assembler  
> **Language:** C++17  
> **Binary:** `assembler_standalone.cpp` → `assembler.exe`

---

## Table of Contents

1. [Problem Statement & Objective](#1-problem-statement--objective)
2. [Project Features](#2-project-features)
3. [Module Diagram, Flowchart & DFD](#3-module-diagram-flowchart--dfd)
4. [Main Modules — Algorithms & Descriptions](#4-main-modules--algorithms--descriptions)
5. [Relevant APS Topics](#5-relevant-aps-topics)
6. [Conclusion](#8-conclusion)

---

## 1. Problem Statement & Objective

### Problem Statement

DNA sequencing machines do not read an entire genome in one pass. Instead, they produce millions of short overlapping fragments called **reads** (typically 75–300 base pairs long). The central challenge of **genome assembly** is to reconstruct the original, complete genomic sequence from these millions of short, noisy, overlapping fragments — a problem analogous to reassembling a shredded book when you have millions of copies of random overlapping excerpts.

This is computationally hard because:
- Reads contain sequencing errors (substitutions, insertions, deletions)
- The genome contains repetitive regions that appear in many places
- The volume of data is enormous (human genome ≈ 3 billion base pairs)
- Both strands of the DNA double helix are sequenced simultaneously

### Objective

Design and implement **HelixForge**, a complete genome assembly pipeline in C++ that:

1. **Ingests** raw sequencing data in FASTQ format
2. **Filters** erroneous k-mers using a Bloom filter (probabilistic data structure)
3. **Constructs** a De Bruijn graph from solid k-mers using a rolling hash
4. **Assembles** the genome using two complementary graph traversal strategies:
   - Hierholzer's Eulerian path algorithm (maximum edge coverage)
   - Dijkstra's shortest path with coverage weights (highest confidence path)
5. **Corrects** sequencing errors in the assembled sequence
6. **Analyses** repeats in the assembly using a Suffix Array and LCP Array
7. **Reports** bioinformatics-standard metrics: GC content, N50, repeat regions

---

## 2. Project Features

| # | Feature | Algorithm / Data Structure | Complexity |
|---|---------|---------------------------|-----------|
| 1 | FASTQ streaming reader | Sequential I/O | O(N) |
| 2 | Reverse complement + canonical k-mers | String manipulation | O(k) per k-mer |
| 3 | Rolling hash (Rabin-Karp) | Polynomial hash with sliding window | O(1) per slide |
| 4 | Solid k-mer filtering | Bloom Filter (dual-filter scheme) | O(1) insert/query |
| 5 | De Bruijn graph construction | Adjacency list + edge frequency map | O(V + E) |
| 6 | Coverage-weighted assembly | Dijkstra's SSSP with min-heap | O((V+E) log V) |
| 7 | Eulerian path assembly | Hierholzer's algorithm | O(E) |
| 8 | Frequency-weighted greedy fallback | Priority-based greedy walk | O(E) |
| 9 | DP error correction | Sliding window majority vote | O(N) |
| 10 | Suffix array construction | Prefix-doubling (Manber & Myers) | O(n log² n) |
| 11 | LCP array construction | Kasai's algorithm | O(n) |
| 12 | Repeat region detection | LCP array scan | O(n) |
| 13 | GC content | Single-pass counter | O(n) |
| 14 | N50 metric | Sort + prefix sum | O(c log c) |
| 15 | FASTA output | Formatted 60-char wrapped output | O(n) |
| 16 | Weighted graph JSON export | Manual JSON writer | O(V + E) |
| 17 | Repeat report (repeats.txt) | Formatted tabular output | O(r) |

**Output Files Generated:**

```
output_dir/
├── genome.fasta       — assembled sequence in FASTA format
├── stats.txt          — full pipeline statistics and timing
├── graph_data.json    — De Bruijn graph (nodes + weighted edges) for visualiser
└── repeats.txt        — repeat region analysis report
```

---

## 3. Module Diagram, Flowchart & DFD

### 3.1 High-Level Module Diagram


<img width="1024" height="1536" alt="image" src="https://github.com/user-attachments/assets/de407279-0af8-41c3-a509-2d8a90996f2b" />


---

### 3.2 Detailed Flowchart


<img width="1024" height="1536" alt="image" src="https://github.com/user-attachments/assets/844f6b3a-e796-40a0-b456-9b1c3896b9de" />


---

### 3.3 Data Flow Diagram (DFD)

#### Level 0 — Context Diagram


<img width="1536" height="1024" alt="image" src="https://github.com/user-attachments/assets/7b679c6f-4421-4b9c-afd4-8dcb90e014ef" />



#### Level 1 — Process Decomposition


<img width="1536" height="1024" alt="image" src="https://github.com/user-attachments/assets/db05352a-0351-4a7c-857e-2757b13209e4" />


#### Level 2 — De Bruijn Graph Data Store


<img width="1536" height="1024" alt="image" src="https://github.com/user-attachments/assets/6c2dbb94-5447-4b73-b70b-77de589c5df2" />


---

## 4. Main Modules — Algorithms & Descriptions

### Module 1 — Reverse Complement (`rev_comp`)

**Purpose:** DNA is double-stranded. A sequencer can read either strand, so the same genomic region may appear as `ATCG` from one read and `CGAT` (its reverse complement) from another. Without canonicalization, these would create duplicate, conflicting graph nodes.

**Algorithm:**
```
rev_comp(s):
    r = ""
    for i from len(s)-1 down to 0:
        complement char: A↔T, C↔G
        append to r
    return r

canonical(kmer):
    return min(kmer, rev_comp(kmer))
```

**Effect:** Every k-mer and its reverse complement now map to the same canonical node, halving graph size and making the graph biologically correct.

---

### Module 2 — Rolling Hash (Rabin-Karp)

**Purpose:** Computing a fresh hash for every k-mer by scanning all k characters is O(k) per k-mer, giving O(Nk) total. Rolling hash reduces each slide to O(1).

**Algorithm:**
```
Polynomial hash:  H(s[i..i+k-1]) = s[i]·B^(k-1) + s[i+1]·B^(k-2) + ... + s[i+k-1]
                  (all mod a Mersenne prime 2^61 - 1)

Slide:  H(s[i+1..i+k]) = (H(s[i..i+k-1]) - s[i]·B^(k-1)) · B + s[i+k]

__uint128_t used for intermediate products to avoid overflow.
```

**Complexity:** O(k) initialisation, O(1) per subsequent slide.

---

### Module 3 — Bloom Filter (Dual-Filter Solid k-mer Detection)

**Purpose:** Remove singleton k-mers caused by sequencing errors. A k-mer seen only once is almost certainly an error, not a real genomic sequence.

**Algorithm:**
```
Two independent Bloom filters: once_BF, twice_BF
Each filter has BITS = 2^26 bits (~8 MB), 3 hash functions.

Insert k-mer hash h:
    IF h ∈ once_BF → insert h into twice_BF
    ELSE           → insert h into once_BF

Query solid: h ∈ twice_BF  (seen at least twice)

Hash function i: bit_pos = ((h XOR salt_i) * knuth_multiplier + offset) mod BITS
```

**False positive rate:** ≈ (1 - e^(-kn/m))^k ≈ 0.02% with k=3, n=10M, m=67M bits.

---

### Module 4 — De Bruijn Graph Construction

**Purpose:** The De Bruijn graph is the core data structure of the assembler. Each (k-1)-mer is a node; each k-mer creates a directed edge from its prefix to its suffix.

**Algorithm:**
```
For each solid k-mer km:
    canonical = min(km, rev_comp(km))
    u = canonical[0..k-2]   (left (k-1)-mer)
    v = canonical[1..k-1]   (right (k-1)-mer)
    adj[u].push_back(v)
    edge_freq["u->v"]++
    outdeg[u]++, indeg[v]++
```

**Key insight:** If a genome has an Eulerian path through its De Bruijn graph, that path spells out the assembled sequence. Euler's theorem guarantees this when every node has balanced in/out degrees.

---

### Module 5 — Dijkstra's Coverage-Weighted Assembly

**Purpose:** Find the highest-confidence assembly path. Edges supported by many reads are more likely to be real. By assigning low cost to high-frequency edges, Dijkstra naturally follows the most-covered route.

**Algorithm:**
```
max_f = max(edge_freq.values())

weight(u→v) = (max_f + 1) - edge_freq[u→v]
    ↑ high coverage = low weight = Dijkstra prefers it

Standard Dijkstra with min-heap (priority_queue, greater<>):
    dist[src] = 0, dist[all others] = ∞
    WHILE heap not empty:
        (d, u) = heap.pop_min()
        IF d > dist[u]: skip  ← lazy deletion of stale entries
        FOR each v in adj[u]:
            nd = d + weight(u→v)
            IF nd < dist[v]:
                dist[v] = nd
                prev[v] = u
                hops[v] = hops[u] + 1
                heap.push((nd, v))

Endpoint: prefer Eulerian sink (indeg - outdeg == 1);
          fallback to node with most hops.

Reconstruct via prev[] back-pointers → sequence.
```

**Complexity:** O((V + E) log V) with binary heap.

---

### Module 6 — Hierholzer's Eulerian Path

**Purpose:** If the De Bruijn graph has an Eulerian path (each node has balanced in/out degrees, at most two exceptions), Hierholzer's algorithm finds a path that traverses every edge exactly once — producing the longest possible assembly.

**Algorithm:**
```
stack.push(start_node)
WHILE stack not empty:
    v = stack.top()
    IF v has unused edges:
        stack.push(next_neighbour[v])
        idx[v]++
    ELSE:
        circuit.append(v)
        stack.pop()
reverse(circuit)

Stitch: seq = circuit[0]
        for i = 1 to len(circuit)-1:
            seq += circuit[i].back()
```

**Complexity:** O(E) — each edge is pushed and popped exactly once.

---

### Module 7 — DP Error Correction

**Purpose:** Sequencing errors often appear as isolated single-base substitutions. A base that disagrees with both its immediate neighbours is almost certainly a substitution error.

**Algorithm:**
```
FOR i = 1 to len(seq)-2:
    IF seq[i-1] == seq[i+1]   ← both neighbours agree
    AND seq[i] ≠ seq[i-1]     ← current base disagrees
        seq[i] = seq[i-1]     ← correct by majority vote
```

**Complexity:** O(N) single pass.  
**DP connection:** Each position's correction decision depends on previously visited positions — this is the optimal substructure property of dynamic programming applied to sequence smoothing.

---

### Module 8 — Suffix Array (Prefix-Doubling)

**Purpose:** The suffix array is the backbone of all string analytics in bioinformatics. It enables O(log n) pattern search, O(n) repeat detection, and is used in production aligners (BWA, Bowtie2).

**Algorithm (Manber & Myers prefix-doubling):**
```
SA  = [0, 1, 2, ..., n-1]      (all suffix start positions)
rank[i] = (unsigned char) s[i]  (initial rank = ASCII value)

FOR gap = 1, 2, 4, 8, ..., n:
    Sort SA using comparator:
        key(i) = (rank[i], rank[i+gap])   (-1 if i+gap out of bounds)
    
    Recompute rank[]:
        rank[SA[0]] = 0
        rank[SA[j]] = rank[SA[j-1]] + (keys differ ? 1 : 0)
    
    IF rank[SA[n-1]] == n-1: BREAK   ← all suffixes uniquely ranked

Sentinel '$' (ASCII 36) appended before build:
    guarantees no two suffixes are identical.
```

**Complexity:** O(log n) doubling rounds × O(n log n) sort = **O(n log² n)** total.

**Why prefix-doubling works:** After round k, rank[] captures the relative order of all length-2^k prefixes. Comparing pairs (rank[i], rank[i+2^k]) gives the order of all length-2^(k+1) prefixes — doubling the resolved length each round.

---

### Module 9 — LCP Array (Kasai's Algorithm)

**Purpose:** The Longest Common Prefix (LCP) array stores the length of the longest common prefix between consecutive suffixes in the suffix array. It is essential for efficient repeat detection.

**Algorithm (Kasai 2001):**
```
Build inverse SA: rank[SA[i]] = i   for all i

h = 0   (current LCP value)
FOR i = 0 to n-1:
    IF rank[i] > 0:
        j = SA[rank[i] - 1]         ← SA-predecessor of suffix i
        WHILE s[i+h] == s[j+h]: h++ ← extend match
        lcp[rank[i]] = h
        IF h > 0: h--               ← key invariant: h drops by ≤ 1
```

**Key invariant:** When we move from suffix i to suffix i+1, the LCP with the SA-predecessor can decrease by at most 1. This means `h` is decremented at most n times total across all iterations → **O(n)** overall.

---

### Module 10 — Repeat Finder (LCP Scan)

**Purpose:** Repetitive DNA regions (transposons, telomeres, microsatellites) confuse assemblers and are biologically significant. The LCP array directly encodes shared prefix lengths, making repeat detection trivial.

**Algorithm:**
```
FOR i = 1 to n-1:
    IF lcp[i] >= min_len:
        Group consecutive entries where lcp[j] >= min_len
        group_lcp = min(lcp[i..j-1])   ← shared prefix of entire group
        Positions: SA[i-1], SA[i], ..., SA[j-1]
        Filter: only positions where pos + group_lcp ≤ seq_len (exclude sentinel)
        Record RepeatRegion{pattern, length, occurrences, positions}

Sort repeat regions by length descending.
```

**Complexity:** O(n) scan of the LCP array.

---

### Module 11 — Greedy Frequency-Weighted Fallback

**Purpose:** If neither Hierholzer nor Dijkstra produces output (highly fragmented graph), a greedy walk provides at least a partial assembly.

**Algorithm:**
```
cur = start_node
res = cur
WHILE cur has unvisited neighbours:
    best_nxt = neighbour v of cur with max edge_freq[cur→v]
               that has not been visited from cur
    res += best_nxt.back()
    cur = best_nxt
```

**Improvement over naive greedy:** Uses edge frequency to prefer the most-supported path at each step, instead of arbitrarily picking the first edge.

---

### Module 12 — Genome Metrics

**GC Content:**
```
gc_content(s) = 100 × |{c ∈ s : c ∈ {G, C}}| / |s|
```
GC content is a fundamental genome characteristic used to identify coding regions, species identity, and DNA stability (G-C pairs have 3 hydrogen bonds vs. A-T's 2).

**N50:**
```
1. Split assembly on 'N' gap characters → contig lengths c[0..m-1]
2. Sort c[] descending
3. Find smallest L such that sum of contigs ≥ L covers ≥ 50% of total assembly
4. N50 = L
```
N50 is the standard assembly quality metric: higher N50 = fewer, longer contigs = better assembly.

---

## 5. Relevant APS Topics

| APS Topic | Where Used in HelixForge | Why It Matters |
|-----------|--------------------------|----------------|
| **Graph Theory — Directed Graphs** | De Bruijn graph construction | Genome assembly reduces to finding a path in a directed graph |
| **Eulerian Paths & Circuits** | Hierholzer's algorithm | An Eulerian path through the De Bruijn graph spells the assembled genome |
| **Shortest Path — SSSP** | Dijkstra's algorithm | Finding the highest-confidence (lowest-cost) path through coverage-weighted graph |
| **Priority Queue / Min-Heap** | Dijkstra implementation | Efficient O(log V) extraction of the minimum-cost node at each step |
| **Hashing** | Rabin-Karp rolling hash | O(1) per k-mer slide instead of O(k) recomputation |
| **Probabilistic Data Structures** | Bloom Filter | Space-efficient approximate set membership for solid k-mer filtering |
| **Divide & Conquer** | Prefix-doubling SA construction | Each round doubles the resolved prefix length; O(log n) rounds total |
| **Sorting** | SA construction (sort inside doubling) | Comparison-based sort drives the O(n log² n) SA algorithm |
| **String Algorithms** | Suffix Array, LCP, Kasai | Fundamental string data structures enabling O(n) repeat detection |
| **Dynamic Programming** | DP error correction; Kasai's LCP invariant | Optimal substructure: each position's state depends on previous |
| **Greedy Algorithms** | Greedy fallback assembler | Locally optimal edge choice (max frequency) at each step |
| **Two-Pointer / Sliding Window** | Rolling hash slide; Kasai's h extension | Amortised O(1) operations using a maintained window |
| **Amortised Analysis** | Kasai's algorithm (h drops ≤ 1 per step) | h increments ≤ n total → O(n) overall despite inner while loop |
| **Space-Time Trade-offs** | Bloom filter (~8 MB vs O(n) exact set) | Accept ≈ 0.02% false positives to save orders of magnitude of memory |
| **Graph Traversal — DFS** | Hierholzer's algorithm (stack-based DFS) | Iterative DFS on the De Bruijn graph finds the Eulerian circuit |
| **Complexity Analysis** | Every module documented | Theoretical vs. measured timing reported in stats.txt |
| **Canonical Forms** | min(kmer, rev_comp(kmer)) | Equivalence classes for bidirectional k-mers; reduces graph size by ~2× |

---

## 6. Conclusion

HelixForge demonstrates that the core algorithms of a real-world genome assembler can be built from first principles using fundamental DAA concepts.

### What We Achieved

1. **Complete end-to-end pipeline** from raw FASTQ reads to assembled FASTA genome, running in under 5 ms on test data.

2. **Two complementary assembly strategies** — Hierholzer (maximum coverage) and Dijkstra (maximum confidence) — are both computed and the better result is used. This dual-strategy approach is novel compared to most textbook implementations.

3. **Bioinformatics-grade string analytics** — the Suffix Array + Kasai LCP implementation correctly identifies all repeat regions, verified against known tandem repeat sequences.

4. **Space efficiency** — the Bloom filter (~8 MB) replaces an exact hash set that would require hundreds of MB for realistic genomes, at a false-positive rate of < 0.02%.

5. **Correctness** — canonical k-mers (reverse complement handling) ensure both DNA strands are treated equivalently, matching the behaviour of production assemblers like SPAdes and Velvet.

### Complexity Summary

| Stage | Algorithm | Complexity |
|-------|-----------|-----------|
| k-mer hashing | Rolling hash + Bloom filter | O(N) |
| Graph construction | Adjacency list | O(V + E) |
| Coverage-weighted assembly | Dijkstra + min-heap | O((V+E) log V) |
| Eulerian assembly | Hierholzer | O(E) |
| Error correction | DP sliding window | O(N) |
| Suffix array | Prefix doubling | O(n log² n) |
| LCP array | Kasai's algorithm | O(n) |
| Repeat detection | LCP scan | O(n) |
| **Overall** | — | **O(N + (V+E) log V + n log² n)** |

### Limitations & Future Work

- **Paired-end reads:** The current pipeline treats all reads as single-end. Paired-end reads provide long-range distance constraints that dramatically improve scaffolding of repeat regions.
- **Contig extraction:** The current implementation outputs one assembled sequence. A full assembler would extract individual unitigs (unbranched paths) as separate contigs.
- **Bubble/tip removal:** Sequencing errors create short dead-end "tips" and two-path "bubbles" in the De Bruijn graph. Removing these before traversal would improve assembly quality.
- **Parallelism:** The hashing stage is embarrassingly parallel — `std::thread` or OpenMP could distribute reads across cores for 4–8× speedup.
- **GFA output:** The Graphical Fragment Assembly format is the industry standard for graph visualisation tools like Bandage.

### Key Takeaways

> The genome assembly problem is a beautiful intersection of graph theory, string algorithms, probabilistic data structures, and dynamic programming — all of which are core DAA topics. HelixForge shows that a competitive implementation requires not just one algorithm but a carefully designed **pipeline** where each module's output feeds the next, with complexity trade-offs made deliberately at every stage.

---

*End of Report — HelixForge Genome Assembler*
