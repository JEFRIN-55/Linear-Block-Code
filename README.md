# Linear Block Code (LBC)

## Aim  
To implement Linear Block Coding (LBC) using a generator matrix to encode messages, calculate Hamming distances, determine error positions, and correct errors in received codewords.

## Tools Required  
- Python  
- NumPy

## Program  
### Linear Block Code Encoding, Syndrome Detection, and Error Correction  
```python
import numpy as np

pb = []  # Parity matrix
Ik = []  # I_K Matrix
p = []
m = []
h = []
h_dis = []
r_code = []
err = []
col = int(input("Enter the Parity bits : "))
row = int(input("Enter the Message bits : "))

# Generator matrix
for i in range(row):
    p = list(map(int, input(f"Enter the row values : {i+1} (Separated by space) : ").split()))
    pb.append(p)
p_mat = np.array(pb, dtype=int)
Ik = np.eye(row, dtype=int)  # Diagonal Matrix
g_mat = np.hstack((p_mat, Ik))  # Generator Matrix

# Codeword length and parity bit length
n, k = g_mat.T.shape

# Possible Message Bits
m = np.array([[1 if (i >> (k - j - 1)) & 1 else 0 for j in range(k)] for i in range(2**k)])

# Codewords and Hamming weights
c = np.mod(np.dot(m, g_mat), 2)

for i, row in enumerate(c):
    h_dis1 = np.sum(row)  # Count number of 1's in the row
    h_dis.append(h_dis1)
h_mat = np.array(h_dis).reshape(1, -1)
d_min = np.min(np.sum(c[1:], axis=1))

# H matrix (Parity-check matrix)
h = p_mat[:, :3]
hp = np.hstack((np.eye(n-k, dtype=int), h.T))
ht = hp.T
zero_row = np.zeros((1, ht.shape[1]), dtype=int)
hp1 = np.vstack((zero_row, ht))

print('**********')
print('The Generator Matrix is: ')
for r in g_mat:
    print(" ".join(map(str, r)))

print('**********')
print(f'Message Bits  Codeword   Hamming Weight')
code_word = np.hstack((m, c, h_mat.T))
for r in range(code_word.shape[0]):
    format_row = " ".join(map(str, code_word[r, :k])) + '\t' + " ".join(map(str, code_word[r, k:n+k])) + '\t' + str(code_word[r, -1])
    print(format_row)

print('**********')
print(f'Minimum Hamming distance : {d_min}')

print('**********')
print(f'Parity Check Matrix')
for r in hp:
    print(" ".join(map(str, r)))

print('**********')
print(f'Parity Check Matrix Transpose')
for r in hp1:
    print(" ".join(map(str, r)))

# Receive codeword
rc = list(map(int, input(f"Enter the error codeword : ").split()))
r_code.append(rc)
r_c = np.array(r_code)

# Syndrome Calculation
e = np.mod(np.dot(r_c, ht), 2)

# Find the Error position
for i in range(n):
    if np.array_equal(e[0], ht[i, :]):
        err = np.eye(n, dtype=int)[i, :]

print(f"The error position is : " + " ".join(map(str, err)))
n1 = hp1.shape[0]
print(f"The size is : {n1}")

print('**********')
print(f'Syndrome Matrix')
for i in range(n1):
    combined_row = np.concatenate((hp1[i, :], np.eye(n1, dtype=int)[i, :]))
    formatted_row = " ".join(map(str, combined_row[:3])) + '\t' + " ".join(map(str, combined_row[k:]))
    print(f'{formatted_row}')

print('**********')
print(f"Syndrome of given received codeword is : " + " ".join(map(str, e[0])))

# Correct the error in the received codeword
add = err + rc
add = np.array(add)
add1 = add % 2
print(f"The correct codeword is : " + " ".join(map(str, add1)))
```

## Output Example  
```
Enter the Parity bits : 3
Enter the Message bits : 3
Enter the row values : 1 (Separated by space) : 1 0 0
Enter the row values : 2 (Separated by space) : 1 0 1
Enter the row values : 3 (Separated by space) : 1 1 0
**********
The Generator Matrix is: 
1 0 0 1 0 0
1 0 1 0 1 0
1 1 0 0 0 1
**********
Message Bits  Codeword     Hamming Weight
0 0 0	      0 0 0 0 0 0	   0
0 0 1	      1 1 0 0 0 1	   3
0 1 0	      1 0 1 0 1 0	   3
0 1 1	      0 1 1 0 1 1	   4
1 0 0	      1 0 0 1 0 0	   2
1 0 1	      0 1 0 1 0 1	   3
1 1 0	      0 0 1 1 1 0	   3
1 1 1	      1 1 1 1 1 1	   6
**********
Minimum Hamming distance : 2
**********
Parity Check Matrix
1 0 0 1 1 1
0 1 0 0 0 1
0 0 1 0 1 0
**********
Parity Check Matrix Transpose
0 0 0
1 0 0
0 1 0
0 0 1
1 0 0
1 0 1
1 1 0
Enter the error codeword : 1 1 0 1 1 1 
The error postion is : 0 0 1 0 0 0
The size is : 7
**********
Syndrome Matrix
0 0 0	  1 0 0 0 0 0 0
1 0 0	  0 1 0 0 0 0 0
0 1 0	  0 0 1 0 0 0 0
0 0 1	  0 0 0 1 0 0 0
1 0 0	  0 0 0 0 1 0 0
1 0 1	  0 0 0 0 0 1 0
1 1 0	  0 0 0 0 0 0 1
**********
Syndeome of given received codeword is : 0 0 1
The correct codeword is : 1 1 1 1 1 1
```
## Manual Calculation
   ![LBC_1(212223060104)](https://github.com/user-attachments/assets/6c587db0-8a26-478e-b476-fbda93afbdf0)
   ![LBC_2(212223060104)](https://github.com/user-attachments/assets/c34c61b4-0bec-45ed-bd94-1d88b8e1c9f2)
   ![LBC_3(212223060104)](https://github.com/user-attachments/assets/dd051afc-c7b9-46b3-bc8c-6b346c8376a4)

## Results  
- The generator matrix was formed using the given parity matrix.  
- All possible message bits and corresponding codewords were generated.  
- The Hamming weights and minimum Hamming distance were computed.  
- The parity-check matrix and its transpose were calculated.  
- The received codeword was verified, and the error position was detected using syndrome decoding.  
- The codeword was successfully corrected.
