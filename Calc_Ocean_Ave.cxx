#include <iostream>
#include<fstream>
#include<numeric>
#include<algorithm>
#include<vector>
#include <sstream>
#define EXIT_FILE_ERROR (1)
#define MAXLINE (10000)
using namespace std;

int main(int argc, char* argv[])
{
  if (argc < 1) {
    cerr << "Usage: " << argv[0] << "requires the species of the absorbing atom !" << endl;
  }
  string prefix1 = "Combo_";
  string prefix2 = "Ave_";
  string ReadFrom_filename = prefix1 + argv[1] + "_all.dat";
  string WriteTo_filename  = prefix2 + argv[1] + "_all.dat";
  vector < vector<double> > Combo_all;
  vector<double> ave_of_row;
  string line;
  ifstream input_file(ReadFrom_filename.c_str());
// make sure file exists and readable 
  if (!input_file.is_open()) {
    cout << "Exititng... unable to open the file" << endl;
    exit(EXIT_FILE_ERROR) ;
  }
  cout << "Reading From file : " << ReadFrom_filename << endl;
// read lines 
  while (getline(input_file, line, '\n')) {
    stringstream ss(line);
    vector<double> rows;
    double value;
    while ( ss >> value ) {rows.push_back(value);}
    Combo_all.push_back(rows);
  }
input_file.close();
/* print out the 2D vecotr to stout
for ( vector< vector<double> >::size_type i = 0; i < Combo_all.size(); i++ ) {
   for ( vector<double>::size_type j = 0; j < Combo_all[i].size(); j++ ) {
      cout << Combo_all[i][j] << ' ';
   }
   cout << endl;
}
*/
//sum up each row to get the combined intensity, then normalize by the number of colums
for ( vector< vector<double> >::size_type i = 0; i < Combo_all.size(); i++ ) {
  double sum_of_elems = std::accumulate(Combo_all[i].begin(), Combo_all[i].end(), 0.0);
  double ave_intensity = sum_of_elems/Combo_all[i].size();
  ave_of_row.push_back(ave_intensity);
}
//write averaged intensity to Ave_O_all.dat
ofstream output_file(WriteTo_filename.c_str());
cout << "Wrting To File: " << WriteTo_filename << endl;
for ( vector<double>::size_type i=0; i<ave_of_row.size(); i++ ) {
  output_file << ave_of_row[i] << endl;
}
output_file.close();
return 0;
}
