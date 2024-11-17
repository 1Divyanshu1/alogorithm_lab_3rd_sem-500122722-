#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

// Function prototypes
void naive_string_match(const char *text, const char *pattern, int *match_count);
void rabin_karp(const char *text, const char *pattern, int *match_count);
void KMPSearch(const char *text, const char *pattern, int *match_count);
void computeLPSArray(const char *pattern, int m, int *lps);

// Naive String Matching Algorithm
void naive_string_match(const char *text, const char *pattern, int *match_count) {
    int n = strlen(text);
    int m = strlen(pattern);

    for (int i = 0; i <= n - m; i++) {
        int j;
        for (j = 0; j < m; j++) {
            if (text[i + j] != pattern[j]) {
                break;
            }
        }
        if (j == m) {
            (*match_count)++;  // Increment match count if pattern is found
        }
    }
}

// Rabin-Karp Algorithm
void rabin_karp(const char *text, const char *pattern, int *match_count) {
    int m = strlen(pattern);
    int n = strlen(text);
    int d = 256;  // Number of characters in the input alphabet
    int q = 101;  // A prime number
    int p = 0;    // Hash value for pattern
    int t = 0;    // Hash value for text
    int h = 1;

    // The value of h would be "pow(d, m-1)%q"
    for (int i = 0; i < m - 1; i++) {
        h = (h * d) % q;
    }

    // Calculate the hash value of pattern and first window of text
    for (int i = 0; i < m; i++) {
        p = (d * p + pattern[i]) % q;
        t = (d * t + text[i]) % q;
    }

    // Slide the pattern over text one by one
    for (int i = 0; i <= n - m; i++) {
        // Check the hash values of current window of text and pattern
        if (p == t) {
            // Check for characters one by one
            int j;
            for (j = 0; j < m; j++) {
                if (text[i + j] != pattern[j]) {
                    break;
                }
            }
            if (j == m) {
                (*match_count)++;  // Increment match count if pattern is found
            }
        }

        // Calculate hash value for next window of text: Remove leading digit, add trailing digit
        if (i < n - m) {
            t = (d * (t - text[i] * h) + text[i + m]) % q;

            // We might get negative value of t, converting it to positive
            if (t < 0) {
                t += q;
            }
        }
    }
}

// Knuth-Morris-Pratt Algorithm
void KMPSearch(const char *text, const char *pattern, int *match_count) {
    int m = strlen(pattern);
    int n = strlen(text);
    int lps[m];  // Longest prefix suffix array

    computeLPSArray(pattern, m, lps);

    int i = 0;  // index for text
    int j = 0;  // index for pattern
    while (i < n) {
        if (pattern[j] == text[i]) {
            i++;
            j++;
        }
        if (j == m) {
            (*match_count)++;  // Increment match count if pattern is found
            j = lps[j - 1];
        } else if (i < n && pattern[j] != text[i]) {
            if (j != 0) {
                j = lps[j - 1];
            } else {
                i++;
            }
        }
    }
}

// Function to compute LPS (Longest Prefix Suffix) array
void computeLPSArray(const char *pattern, int m, int *lps) {
    int len = 0;  // Length of the previous longest prefix suffix
    lps[0] = 0;   // lps[0] is always 0
    int i = 1;

    while (i < m) {
        if (pattern[i] == pattern[len]) {
            len++;
            lps[i] = len;
            i++;
        } else {
            if (len != 0) {
                len = lps[len - 1];
            } else {
                lps[i] = 0;
                i++;
            }
        }
    }
}

// Function to measure time taken by an algorithm
double measure_time(void (*func)(const char*, const char*, int*), const char *text, const char *pattern, int *match_count, int repeats) {
    clock_t start, end;
    start = clock();
    for (int i = 0; i < repeats; i++) {
        func(text, pattern, match_count);
    }
    end = clock();
    return ((double)(end - start)) / CLOCKS_PER_SEC; // Return time in seconds
}

// Main function to test the algorithms
int main() {
    const char *patterns[5] = {
        "ABAB",
        "ABCABCDABAB",
        "AABCAAADA",
        "ABABABABABABABABAB",
        "ABCDABCDABCDABCDABCDABCD"
    };

    const char *texts[5] = {
        "ABABDABACDABABCABABABABDABACDABABCABABABABDABACDABABCABABABABABDABACDABABABABDABACDABABABABABDABACDABABCABABABAB",
        "AABABAAABABAAABABAAABABAAABABAAABABAAABABAAABABABABABABABAAABABABABABABABABABABABABABABABABABABABABABAB",
        "AABAAABAAABAAABAAABAAABAAABAAABAAABAAABAAABAAABAAABAAABAAABAAABAAABAAABAAABAAABAAABAAABAAABAAABAAABAA",
        "ABABABABABABABABABABABABABABABABABABABABABABABABABABABABABABABABABABABABABABABABABABABABABABABABABABABABAB",
        "ABCDABCDABCDABCDABCDABCDABCDABCDABCDABCDABCDABCDABCDABCDABCDABCDABCDABCDABCDABCDABCDABCDABCDABCDABCDABCD"
    };

    int repeats = 1000;  // Number of times to repeat each algorithm

    // Loop through each input size
    for (int i = 0; i < 5; i++) {
        printf("Input %d:\n", i + 1);
        int match_count_naive = 0, match_count_rabin = 0, match_count_kmp = 0;

        // Measure and print time for Naive String Matching
        double time_naive = measure_time(naive_string_match, texts[i], patterns[i], &match_count_naive, repeats);
        printf("Naive String Matching:\nTime taken: %f seconds\nMatches found: %d\n\n", time_naive, match_count_naive);

        // Measure and print time for Rabin-Karp Algorithm
        double time_rabin = measure_time(rabin_karp, texts[i], patterns[i], &match_count_rabin, repeats);
        printf("Rabin-Karp Algorithm:\nTime taken: %f seconds\nMatches found: %d\n\n", time_rabin, match_count_rabin);

        // Measure and print time for Knuth-Morris-Pratt Algorithm
        double time_kmp = measure_time(KMPSearch, texts[i], patterns[i], &match_count_kmp, repeats);
        printf("Knuth-Morris-Pratt Algorithm:\nTime taken: %f seconds\nMatches found: %d\n\n", time_kmp, match_count_kmp);

        // Final Comparison for the input
        printf("Performance Comparison for Input %d:\n", i + 1);
        printf("Naive: %f seconds, Matches: %d\n", time_naive, match_count_naive);
        printf("Rabin-Karp: %f seconds, Matches: %d\n", time_rabin, match_count_rabin);
        printf("KMP: %f seconds, Matches: %d\n\n", time_kmp, match_count_kmp);
    }

    return 0;
}

