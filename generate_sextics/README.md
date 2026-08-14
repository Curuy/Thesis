
Generating $S_3$ sextic polynomials in a given discriminant range `[Y, X]`. 

It uses **CUBIC** package (v1.4) developed by Karim Belabas.
---

## 1. Prerequisites


*   **GCC** (with OpenMP support for multithreading)
*   **PARI/GP** library 
*   **Make**

*Note: Ensure your PARI library is installed in a standard location (e.g., `/usr` or `/usr/local`). If it is installed elsewhere, update the `PARIHOME` variable at the top of the `Makefile`.*

---

## 2. Compilation

To compile the project, open your terminal in the project root and run:

```bash
make

### Builds
*   **`make noprint`**: (Recommended for large ranges). Compiles an optimized version that does not print every matched field to the terminal. This is crucial for heavily multithreaded runs to avoid terminal I/O bottlenecks and segmentation faults.
*   **`make debug`**: Compiles with debug flags enabled.
*   **`make profile`**: Compiles with profiling flags (`-pg`) for performance analysis.
*   **`make clean`**: Safely removes all compiled object files, test binaries, and build directories.

---

## 3. Usage

The compiled executable is named `sextic` and accepts mathematical expressions (like `10^10`) natively using PARI's string parser.

**Syntax:**
```bash
./sextic <Y> <X> [num_threads]
```
*   `Y`: Starting bound.
*   `X`: Ending bound 
*   `num_threads` *(Optional)*: Number of threads. Default is use all threads 
**Examples:**
```bash
# Search from 0 to 10^10 all available threads
./sextic 0 10^10

# Search from 10^12 to 10^14 using 4 threads
./sextic 10^12 10^14 4
```

---

## 4. Running the Tests

Has some tests which compare results from 
https://arxiv.org/pdf/1208.2170

**To run the entire test suite automatically:**
```bash
./run_tests.sh
```

AI was used to generate this README

**Author of original algorithm:**
Karim Belabas, IMB (UMR 5251)
Univ. Bordeaux, 351 cours de la Liberation, F-33405 Talence (France)
http://www.math.u-bordeaux.fr/~kbelabas/

**References:**

> Belabas, Karim. "A fast algorithm to compute cubic fields." *Mathematics of Computation* 66.219 (1997): 1213-1237.
> MR1415795 (97m:11159)

> Belabas, Karim. "On quadratic fields with large 3-rank." *Mathematics of Computation* 73.248 (2004): 2061-2074.
> MR2059751 (2005c:11132)
