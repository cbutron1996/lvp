This program runs the Landau-Vishkin 89 (LV-89) algorithm using SDSL's Compressed Suffix Trees

LV 89 is an approximate pattern matching algorithm that finds all occurrences of a PATTERN in a TEXT up to k errors. It accepts text files for fill in the TEXT, strings to fill in the PATTERN, and an integer to fill in the k errors.

The program comprises of only one file LV.cpp. The other files are just input text files.

How to Compile: g++ -std=c++11 -O3 -DNDEBUG -I ~/include -L ~/lib LV.cpp -o ./LV.out -lsdsl -ldivsufsort -ldivsufsort64
How to Run: ./LV.out [integer k errors] [string text file] [string pattern]

Test Results:
Input: ./LV.out 2 text4KB aaa
- text4KB - a 5 KB text file repeating the word "baananaaan"
- pattern - "aaa"
- k = 2
Output:
-	Num Occurrences: 9,117
-	Run Time: 0.140 seconds (140 milliseconds)

Input: ./LV.out 2 text1MB aaa
- text1MB - a 1.1 MB text file repeating the word "baananaaan"
- pattern - "aaa"
- k = 2
Output:
-	Num Occurrences: 2,211,837
-	Peak Memory Usage: 10 MB
-	Run Time: 15.869 seconds

Input: ./LV.out 2 text18MB aaa
- text18MB - a 17.7 MB tex file repeating the word "baananaaan"
- pattern - "aaa"
- k = 2
Output:
-	Num Occurrences: 35,389,437
-	Peak Memory Usage: 147 MB
-	Run Time: 250 seconds

Input: ./LV.out 3 simpleText ABCAAB
- simpleText - an 19 byte text file with the text "ABABCAABCAABCAABAA"
- pattern - "ABCAAB"
- k = 3
Output:
-	Num Occurrences: 40
-	Run Time: 0.087 seconds (87 milliseconds)
