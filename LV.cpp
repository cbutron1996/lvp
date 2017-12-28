#include <sdsl/suffix_trees.hpp>
#include <iostream>
#include <fstream>
#include <string>
using namespace sdsl;
using namespace std;

using namespace std::chrono;
using timer = std::chrono::high_resolution_clock;

class MatrixL {
private:
	int k;
	int n;
	int m;
	int **matrix; //Matrix
	string str;
	cst_sada<> cst; //CST
public:
	string concat(string s1, string s2) {
		return s1 + "#" + s2;
	}
	MatrixL(int maxErrors, string text, string pattern) {

	    // ifstream infile(file);
	    // string text( (istreambuf_iterator<char>(infile)), (istreambuf_iterator<char>()) );
	    // //string text = "bananaaan";
	    // cout << "Finish reading data" << endl;

	    //ofstream outfile("outputText.txt", ofstream::out);
	    //outfile << text;
	    //cout << "Finish writing data" << endl;

		k=maxErrors;
		n=text.length();
    m=pattern.length();

		str = concat(text, pattern);

		construct_im(cst, str, 1); // new CST version

		cout << "Finish constructing cst" << endl;

		//the LCP array is built in this
		//SuffixArray constructor

        	//matrix = new int [n-m+k+1][k+1] //from java
        	int numRows = n-m+k+1;
        	int numCols = k+1;
    		matrix = new int *[numRows];
    		for (int i = 0; i < numRows; ++i)
            		matrix[i] =new int[numCols];

        	//init to all 0's YL
		for(int i = 0; i < numRows; i ++){
        		for(int j = 0; j < numCols; j++){
                		matrix[i][j] = 0;
                }
        }
        cout << "Finished step 1" << endl;
        // computeMatrixL();

				initStepTwo();
    		cout << "Finished step 2" << endl;
    		initStepThree();
    		cout << "Finished step 3" << endl;
    		fillMatrix();
    		cout << "Filled matrix" << endl;
	}

	~MatrixL() {
		int numRows = n-m+k+1;
		for (int i = 0; i < numRows; i++) {
        		delete [] matrix[i];
    		}
    		delete [] matrix;
	}

	// void computeMatrixL() {
  //   		initStepTwo();
  //   		cout << "Finished step 2" << endl;
  //   		initStepThree();
  //   		cout << "Finished step 3" << endl;
  //   		fillMatrix();
  //   		cout << "Filled matrix" << endl;
  //   		printMatrix();
  //   		cout << "Printed matrix" << endl;
  //   	}

	void initStepTwo(){
    		int d, i, j=k;
    		int row,col;

    		for (d=-(j), i=1; d < -1 &&  i>=1; d++, i++){
          	  	row = transform(d);
            		col = j-i;
        		//matrix[transform(d)][j-i]=j-i;
        		matrix[row][col]=j-i;
    		}
    	}

	void initStepThree () {
    		int d, i, j=k;
    		int row,col;

    		for (d=-(j), i=1; d < 0 &&  i>=1; d++, i++){
       			row = transform(d);
       			col = j-i+1;
        		// matrix[transform(d)][j-i+1]=j-i+1;
        		matrix[transform(d)][j-i+1]=j-i+1;
    		}
    	}

	int transform (int x) {
    		return x+k;
    	}

	void fillMatrix () {
		int numOccurrences = 0;

    int e, d, row, upperLeft, immediateLeft, lowerLeft, indexOfOccurence=-1;
    bool isBottomRow=false;
    upperLeft=immediateLeft=lowerLeft=-1;

    for (e=0; e<=k; e++){ //iterate over each column
    	for (d= (e*-1)+1; d <= n-m; d++){ //iterate over select rows

 		   //#cout<< "\nlogical d: " << d << " physical d: " << transform(d) << " e: " << e << endl;

 		   //all elements in the first physical column have their left values initialized to -1,
 		   //as specified in the first initialization step of the Landau-Vishkin paper
 		   	if (e==0){
 				upperLeft=immediateLeft=lowerLeft=-1;
 				//#cout << "Immediate left (initialized): " << immediateLeft << endl;
 				//#cout << "Lower left (initialized): " << lowerLeft << endl;
 				//#cout << "Upper left (initialized): " << upperLeft << endl;

 			}

     		   	//all elements not in the first physical column must actually look for their left values
 			else{

				immediateLeft=matrix[transform(d)][e-1];
				//#cout << "Immediate left: " << matrix[transform(d)][e-1] << endl;

				upperLeft=matrix[transform(d-1)][e-1];
				//#cout << "Upper left: " << matrix[transform(d-1)][e-1] << endl;


 				if(transform(d)<transform(n-m)){//if the element is NOT in the bottom row

 					lowerLeft=matrix[transform(d+1)][e-1]; //its lower left value can be read
     				//#cout << "Lower left: " << matrix[transform(d+1)][e-1] << endl;
 				}

     				//if the element is in the bottom row, its lower left element obviously cannot be read
 				else{

 					//#cout << "Cannot read from spot to bottom left at: " << transform(d+1) << " ";
 					//#cout << e-1 << endl;
 					isBottomRow=true;
 				}

 			}//end else


     		   	//if the element is on the bottom row, then row is computed as the max of only two left values
 			if (isBottomRow==true){//changed to == YL
 				row = max(upperLeft+1, immediateLeft+1);
 				//#cout << "max: " << row << endl;
 				isBottomRow=false;
 			}

 			//otherwise, row is computed by looking at three left values
 			else{
				 row =  max(max(immediateLeft+1,lowerLeft),upperLeft+1);
				 //#cout << "max: " << row << endl;
 			}

 			row = min(row, m);

 			//#cout << "m: " << m << endl;
 			//#cout << "row: " << row << endl;

 			//int lcp= suffixArray->calculateLCP (str, row+d, row+n+1, n, m);
 			// From RMQ to CST LCA
			int rank1 = cst.csa.isa[row+d]; // first node
			int rank2 = cst.csa.isa[row+n+1]; // second node
      auto ancestor = cst.node(rank1,rank2); // works in cst_sada
			// auto ancestor = cst.lca(cst.node(first,first), cst.node(second,second)); // works in others
			int lcp = cst.depth(ancestor); // CST Depth of LCA calls = LCP

			//row+n+1 is the index in str that corresponds to the index of row in pattern
 			//#cout << "lcp: " << lcp << endl;

 			matrix[transform(d)][e]= row + lcp;
 			//#cout << "matrix element: " << matrix[transform(d)][e] << endl;

 			if ( (matrix[transform(d)][e] == m) && (d+m<=n) ){
 				indexOfOccurence=d+m-1;
 				// cout << "There is an occurrence ending at index " << indexOfOccurence << " of the text when e=" << e << endl;
				numOccurrences++;
			}



   		}//end for

   	}//end for

  	//if (indexOfOccurence==-1)
  	//cout << "No occurrences" << endl;
		cout << "Num Occurrences: " << numOccurrences << endl;
	}//end fillMatrixL

	//prints matrix
	void printMatrix() {
        	int numRows;
        	int tempInt = n-m+k+1;
        	numRows = tempInt;
        	int numCols = k+1;
        	for(int i = 0; i < numRows; i ++){
            		for(int j = 0; j < numCols; j++){
                		cout << matrix[i][j] << " ";
            		}
            		cout << endl;
        	}
	}
};


int main(int argc, char* argv[]) {
    if (argc < 4) {
        cout << "Usage: " << argv[0] << " <maxErrors> <file> <pattern>" << endl;
        return 1;
    }

		ifstream infile(argv[2]);
		string text( (istreambuf_iterator<char>(infile)), (istreambuf_iterator<char>()) );
		//string text = "bananaaan";
		cout << "Finish reading data" << endl;

    auto start = timer::now();
	MatrixL matrixL = MatrixL(strtol(argv[1], NULL, 10), text, argv[3]);
	auto stop = timer::now();
	cout << "Program Run: Time: " << duration_cast<seconds>(stop-start).count() << endl;

	std::cout << "peak usage = " << memory_monitor::peak() / (1024) << " KB" << std::endl;
	std::ofstream cstofs("cst-construction.html");
    memory_monitor::write_memory_log<HTML_FORMAT>(cstofs);
    cstofs.close();

	return 0;
}
