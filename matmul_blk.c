void matmul_blk(int m, int n, int k, double **A, double **B, double **C, int bs){
	// Loop iterators
	int r, s, t;
	int rr, ss, tt;
	int p = m/bs; //number of rows of submatrices in A and C
 	int q = n/bs; //number of columns of submatrices in B and C
	int u = k/bs; //number of columns in A and rows in B, i.e. "shared dimension" of A and B

	for(r=0; r < m; r++){
		for(s=0; s < n; s++){
			C[r][s] = 0;
		}
	}

	// First loop through all submatrices of C
	for(rr=0; rr<p; rr++){ //iterate over rows of submatrices
		for(ss=0; ss<q; ss++){ //iterate over columns of submatrices
			for(tt=0; tt<u; tt++){ //compute matrixproduct of submatrices of row rr and column ss and do this tt times
				// Compute submatrix product as usual on all submatrices in the tt-tt dimensions
				for(r=0; r<bs; r++){
					for(s=0; s<bs; s++){
						for(t=0; t<bs; t++){
							C[rr*bs+r][ss*bs+s] = C[rr*bs+r][ss*bs+s] + A[rr*bs+r][tt*bs+t]*B[tt*bs+t][ss*bs+s];
						} 
					}
				}
			}
		}
	}
}
