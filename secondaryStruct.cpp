#include <iostream>
#include <vector>
#include <set>

using namespace std;

/**
 * Calculates the secondary structure of an RNA string using dynamic programming.
 *
 * @param RNA the input RNA string
 *
 * @return a set of pairs representing the start and end indices of the secondary structure
 *
 * @throws None
 */
set<pair<int, int>> secondaryStruct(string &RNA)
{
	time_t start, end;
	time(&start);
	// Initialise the dynamic programming vector
	int n = RNA.size();
	vector<vector<int>> opt(n, vector<int>(n, 0));

	// Initialise the base cases
	if (n < 4)
		return set<pair<int, int>>();

	// Calculate the maximum number of base pairs in the secondary structure
	for (int k = 4; k <= n; k++)
	{
		for (int i = 0; i <= n - k; i++)
		{
			int j = i + k - 1;
			// Check if the base pairs are complementary
			if (RNA[i] == 'A' && RNA[j] == 'U' || RNA[i] == 'U' && RNA[j] == 'A' || RNA[i] == 'C' && RNA[j] == 'G' || RNA[i] == 'G' && RNA[j] == 'C')
				opt[i][j] = opt[i + 1][j - 1] + 1;

			// Check if the base pairs are not complementary
			for (int k = i; k < j; k++)
				opt[i][j] = max(opt[i][j], opt[i][k] + opt[k + 1][j]);
		}
	}

	set<pair<int, int>> structure;
	set<pair<int, int>> toAdd;
	for (int i = 0; i < n; i++)
	{
		for (int j = i + 4; j < n; j++)
		{
			if (((RNA[i] == 'A' && RNA[j] == 'U') || (RNA[i] == 'U' && RNA[j] == 'A') || (RNA[i] == 'C' && RNA[j] == 'G') || (RNA[i] == 'G' && RNA[j] == 'C')) && (opt[i][j] == opt[i + 1][j - 1] + 1))
			{
				bool valid = true;
				for (auto &pair : toAdd)
				{
					int k = pair.first, l = pair.second;
					if ((i < k && k < j && j < l) || (i == k && j < l) || (i < k && k < j && j == l) || (i < k && k < j && j > l))
					{
						valid = false;
						break;
					}
				}
				if (valid)
					toAdd.insert({i, j});
			}
		}
	}

	// cout << "toAdd Structure:";
	// for (const auto &pair : toAdd)
	// 	cout << "(" << pair.first << "," << pair.second << ") ";

	for (auto &pair : toAdd)
	{
		bool flag = true;
		for (auto &existing : structure)
		{
			int k = existing.first, l = existing.second;
			if ((pair.first > k && pair.second == l) || (pair.first > k && l > pair.first && pair.second > l))
			{
				flag = false;
				break;
			}
			if ((pair.first == k && pair.second > l) || (pair.first < k && pair.second == l))
			{
				structure.erase(existing);
				break;
			}
		}
		if (flag)
			structure.insert(pair);
	}

	time(&end);
	// Calculating total time taken by the program.
	long double time_taken = end - start;
	// cout << "Time taken by program is : " << fixed << time_taken;
	// cout << " sec " << endl;
	return structure;
}

/**
 * Reads an RNA string from the standard input, calculates its secondary structure using dynamic programming,
 * and prints the start and end indices of the secondary structure to the standard output.
 *
 * @return 0 if the function executes successfully
 *
 * @throws None
 */
int main()
{
	// Read an RNA string from the standard input
	string RNA;
	cin >> RNA;

	// Calculate the secondary structure using dynamic programming
	set<pair<int, int>> structure = secondaryStruct(RNA);

	// Print the start and end indices of the secondary structure
	cout << "Secondary Structure: ";
	for (const auto &pair : structure)
		cout << "(" << pair.first << "," << pair.second << ") ";

	cout << endl;
	return 0;
}
