#include <iostream>
#include <sstream>
#include <string>
#include <vector>
#include <unordered_map>
#include <algorithm>

int main() {
    std::string line;
    std::unordered_map<std::string,std::vector<int> > quals;
    std::unordered_map<std::string,int > quals_threshold;
    std::vector<std::string> samples_genotypes;
    std::vector<int> samples_genotypes_quality;

    samples_genotypes.reserve(5000);
    samples_genotypes_quality.reserve(5000);

    // Process each line from stdin
    while (std::getline(std::cin, line)) {
        // Check if the line is a header line (starts with '#')
        if (line[0] == '#') {
            std::cout << line <<"\n"; // Print header as it is
        } else {
            // Split VCF line for computations
            std::istringstream iss(line);

            for(auto q:quals)
            {
                q.second.clear();
            }
            samples_genotypes.clear();
            samples_genotypes_quality.clear();

            std::string token;
            int field_count=0;
            while (std::getline(iss, token, '\t')) { // Split by tab character
                if(field_count < 9) {
                    std::cout << token << "\t";
                }
                else{
                    unsigned firstColon = token.find(':');
                    unsigned secondColon = token.find(':',firstColon+1);
                    std::string genotype=token.substr(0,firstColon);
                    int quality =  atoi(token.substr(firstColon+1,secondColon-firstColon).c_str());

                    samples_genotypes.push_back(genotype);
                    samples_genotypes_quality.push_back(quality);

                    auto quals_iter = quals.find(genotype);
                    if(quals_iter == quals.end())
                    {
                        quals[genotype] = std::vector<int>();
                        quals_iter = quals.find(genotype);
                        quals_iter->second.reserve(5000);
                    }
                    quals_iter->second.push_back(quality);

                }
                field_count++;
            }
            for(auto q:quals)
            {
                sort(q.second.begin(),q.second.end());
                quals_threshold[q.first] = q.second[q.second.size()/2];
            }

            for(unsigned i=0;i< samples_genotypes.size();i++)
            {
                if(samples_genotypes_quality[i] < quals_threshold[samples_genotypes[i]])
                    std::cout << "./.:"<<samples_genotypes_quality[i];
                else
                    std::cout <<samples_genotypes[i] << ":"<<samples_genotypes_quality[i];
                if(i!= samples_genotypes.size()-1)
                    std::cout << "\t";

            }



            std::cout << "\n";

            // Add your computations here
        }
    }

    return 0;
}
