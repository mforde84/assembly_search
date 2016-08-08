#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <thread>
#include <math.h>
#define o std::cout
#define s std::string 
#define v std::vector
#define ss search_space
#define i std::ifstream
#define t std::thread

//object to hold annotation specific consensus sequence
struct search_space{
	
	s annotation, space;

};

v<ss> load_consensus(auto *file){

	//read line delimited text file holding consensus sequence
	//push to vector of search space objects
	//each new object defined by > character indicating a chr or patch annotation.
	int annotation_count = -1;
	s line;
	v<ss> push;
	ss holder;
	i input(file);
	while (getline(input,line)){
		if(line.at(0) == '>'){ 
			//if hit annotation line, push empty struct, and assign annotation to object
			annotation_count++;
			push.push_back(holder);
			push.at(annotation_count).annotation = line;
		}else 
			//else concatenate line to string holding consensus, overflow at ~4.2B char
			push.at(annotation_count).space += line;
	}
	input.close();
	return(push);

}

v<s> load_query(auto *file){

	//read line delimited file with search queries
	//assign vector of strings
	s line;
	v<s> push;
	i input(file);
	while (getline(input,line))
		push.push_back(line);
	input.close();
	return(push);

}

void print_hits(auto &assembly, auto &query, auto *file){

	//greedy search 
	//iterate through each ss object for each query string 
	for(query_token : query){
		for(assembly_token : assembly){
			//find first instance of query
			std::size_t found = assembly_token.space.find(query_token);
			if(found != s ::npos){
				//hits to stout
				o << file << "\t" << assembly_token.annotation << "\t" << query_token << "\t" << found << "\n"; 
				bool hit = 1;
				while (hit){
					//find additional instances of query in ss consensus object offset by 1 
					found = assembly_token.space.find(query_token,found+1);
					if(found != s ::npos)
						//hits to stout
						o << file << "\t" << assembly_token.annotation << "\t" << query_token << "\t" << found << "\n"; 
					else
						hit = 0;
				}
			}
		}
	}			

}

void runtime_threads(char *file, v<s> &query){

	//calls functions to load assembly to memory, and output query hits to stout
	v<ss> assembly = load_consensus(file);
	print_hits(assembly, query, file);

}

void define_threads(char *argv[], int n_samples, int n_threads, auto &query){
	
	//instantiates then joins threads
	v <t> thread_holder;
	if(n_samples <= n_threads){
		//if number of samples is less than or equal to num of threads, then just
		//run a single sample per thread
		for (int n = 0; n < n_samples; n++)
			thread_holder.push_back(t(runtime_threads,argv[3+n],std::ref(query)));
		for (threads : thread_holder)
			threads.join();
	}else{
		//else iterate through subsets of threads
		double n_iterations = floor((n_samples / n_threads) + 0.5);
		int n_per_thread = n_samples / n_threads, n_leftover = n_samples % n_threads, mod = n_per_thread;
		for (int n = 0; n < n_iterations; n++){
			if (n == (n_iterations - 1))
				mod += n_leftover;
			for (int x = 0; x < mod; x++){
				int ptr = (n*n_per_thread)+x;
				thread_holder.push_back(t(runtime_threads,argv[ptr],std::ref(query)));
			}
			for (threads : thread_holder)
				threads.join();
		}
	}
	
}

int main (int argc, char* argv[]) { 
	
	//usage: $ assembly-search n-threads dictionary.txt consensus-assembly.txt > output.txt 
	int n_samples = argc - 3, n_threads = atoi(argv[1]);
	v<s> query = load_query(argv[2]);
	define_threads(argv, n_samples, n_threads, query);
	
}

