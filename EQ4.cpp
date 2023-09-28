#include <iostream>
#include <vector>

// Description:
// The main function implements Algorithm T (lexicographic combinations) from "The Art of Computer Programming"
// by Donald Knuth [found in Volume 4A on page 359] to generate all k-combinations from a set of J elements {1,...,J}. 
// The function then prints each combination.
// Accepts:
// int k: Represents the size of each combination.
// int J: Represents the size of the set from which to pull elements.
// Returns:
// Prints each combination
// int 0: If the program executed successfully.
// int 1: If invalid input is provided.

int main() {
    // Declarations of k, J to store user input
    int k, J;
    std::cout << "Enter the value of k: ";
    std::cin >> k;
    std::cout << "Enter the value of J: ";
    std::cin >> J;

    // Check for invalid input
    if(k >= J || k <= 0 || J <= 0) {
        std::cerr << "Invalid input. Please ensure 0 < k < J." << std::endl;
        return 1; // Return 1 indicating an error occurred due to invalid input
    }

    // Initialize vector c with (k + 3) elements.
    std::vector<int> c(k + 3);
    
    // Step T1: Initialize c to hold the first combination {1, 2, ..., k}
    for(int j = 1; j <= k; ++j)
        c[j] = j; 

    // Setup sentinel values. c[k + 1] acts as a boundary sentinel.
    c[k + 1] = J + 1;
    // c[k + 2] holds the value 0 and is not used in generating combinations but is required as per Algorithm T
    c[k + 2] = 0;
    
    // j is the smallest index such that c[j + 1] > j, initialized to k.
    int j = k; 

    while (true) {
        // Step T2: Visit. Print the current combination.
        for(int i = k; i >= 1; --i)
            std::cout << c[i] << ' ';
        std::cout << '\n';

        // Exit the loop if j == 0, meaning all combinations have been visited
        if(j == 0) break;
        
        // Store the value of j in x before going to Step T6.
        int x = j;

        // Step T3: Easy Case? Increment c[1] if it can be incremented without exceeding c[2].
        if(c[1] + 1 < c[2]) {
            ++c[1];
            continue;
        } else { // If not, proceed to find the rightmost position to increment.
            j = 2;
        }

        // Step T4: Find j. This step finds the rightmost index that can be incremented.
        while (true) {
            // Update c[j - 1] before going to the next index.
            c[j - 1] = j - 1;
            x = c[j] + 1;
            // Exit this loop if we found the rightmost position to increment.
            if(x < c[j + 1]) break;
            // Move to the next index to the right.
            ++j;
            // If j > k, we have visited all combinations, and we terminate the algorithm.
            if(j > k) return 0;
        }

        // Step T6: Decrease c_j. Update c[j] and move back to the previous index.
        c[j] = x;
        --j;
    }

    // Return 0 indicating that the program executed successfully.
    return 0;
}
