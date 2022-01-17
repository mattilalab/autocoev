#include<iostream>
#include<fstream>
#include"BCFasta.h"
#include"file_manip.h"
#include"caps.h"
#include"create.h"
#include<getopt.h>
#include<math.h>
#include<gsl/gsl_randist.h>
#include<gsl/gsl_rng.h>
#include<gsl/gsl_cdf.h>
#include<gsl/gsl_matrix.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_statistics.h>
#include<sys/time.h>
#include<iomanip>
#include <bits/stdc++.h> 



#include <Seq/SequenceApplicationTools.h>
#include <Seq/SiteTools.h>
#include <Seq/Sequence.h>
#include <Seq/SequenceTools.h>
#include <Seq/ioseq>
#include <Seq/alphabets>

// From Utils:
#include <Utils/ApplicationTools.h>
#include <Utils/TextTools.h>
#include <Utils/KeyvalTools.h>


// From PhylLib:
#include <Phyl/Tree.h>
#include <Phyl/Node.h>
#include <Phyl/DRHomogeneousTreeLikelihood.h>
#include <Phyl/DRNonHomogeneousTreeLikelihood.h>
#include <Phyl/PatternTools.h>
#include <Phyl/PhylogeneticsApplicationTools.h>
#include <Phyl/MarginalAncestralStateReconstruction.h>
#include <Phyl/OptimizationTools.h>
#include <Phyl/RASTools.h>
#include <Phyl/TreeLikelihoodTools.h>
#include <Phyl/TreeTemplateTools.h>
#include <Phyl/DRTreeLikelihoodTools.h>
#include <Phyl/MarkovModulatedSubstitutionModel.h>
#include <Phyl/SubstitutionModelSet.h>
#include <Phyl/SubstitutionModelSetTools.h>
#include <Phyl/Newick.h>

// From NumCalc:
#include <NumCalc/DiscreteDistribution.h>
#include <NumCalc/ConstantDistribution.h>
#include <NumCalc/DataTable.h>
#include <NumCalc/MatrixTools.h>
#include <NumCalc/VectorTools.h>
#include <NumCalc/AutoParameter.h>

// From Utils:
#include <Utils/BppApplication.h>
#include <Utils/ApplicationTools.h>
#include <Utils/FileTools.h>
#include <Utils/TextTools.h>



//static int statl=0;
const gsl_rng_type * T;
gsl_rng *r;

vector<double> totaltempnew;
double alphathresh = 0;
int main(int argc, char *argv[]){



	string pwd, temp_folder, mat_file, aminos, folder, forced, tree_file, mystring, structure, output, treefile, pdb_folder;
	double threshold, thresholdR, gapthresh=0.5, threshval=0.01, bootcut=0.75;
	int num_randoms, three, time_corr, analysis=0, opt=0, seed, variable=0, tree_in=0, pdb_true=0, converge=0;
	struct timeval start;

	print_splash("stdout");

	if(argc<3){
		fprintf(stderr, "usage: [-F alignmentfolder] options([-S Structure folder]\n[--intra/--inter analysis type][-a alpha value][-r random samples][-N Newick formatted tree]\n");
		exit(-1);
	}

	Pwd(pwd);

	threshold=0.01;num_randoms = 100;thresholdR = 0.01;analysis=0;three=1;time_corr=0;

	/*	struct globalArgs_t {
		int noIndex;     
		char *langCode; 
		const char *outFileName;
		FILE *outFile;
		int verbosity;         
		char **inputFiles;    
		int numInputFiles;   
		int randomized;     
		} globalArgs;
		*/
	string refname;
	static const struct option longOpts[] = {
		{ "inter", no_argument, NULL, 'i' },
		{ "intra", no_argument, NULL, 'l' },
		{ NULL, no_argument, NULL, 0 }
	};

	bool treefol=false;
	int *longIndex=NULL;

	while ((opt = getopt_long(argc, argv, "H:b:cF:o:r:ilg:a:G:T:vS:R:", longOpts, longIndex)) != -1) {
		switch(opt) {

			case 'T':
				treefol = true;
				treefile = optarg;
				break;
			case 'v':
				variable = 1;
				break;
			case 'i':
				analysis = 1;
				break;
			case 'l':
				analysis = 0;
				break;
			case 'a':
				threshval = atof(optarg);
				break;
			case 'S':
				pdb_true=1;
				pdb_folder = optarg;
				break;
			case 'F':
				if(optarg[0]=='~'||optarg[0]=='/'){
					mystring = optarg;
					break;
				}else{
					char temp[1000];
					getcwd(temp, 1000);
					mystring += temp;
					mystring  += "/";
					mystring += optarg;
					mystring += "/";
					break;
				}
			case 'R':
				thresholdR=atof(optarg);
				break;
			case 'b':
				bootcut = atof(optarg);
				break;
			case 'r':
				num_randoms = atoi(optarg);	
				break;
			case 'g':
				gapthresh=atof(optarg);
				break;
			case 'H':
				refname = optarg;
				break;
			case 'c':
				converge = 1;
				break;
			case ':':
				fprintf(stderr, "Error: Unknown option passed in: %c\n", optopt); /* optarg defined in getopt.h */
				fprintf(stderr, "usage: [-F alignmentfolder] options([-S Structure folder]\n[--intra/--inter analysis type][-a threshold alpha value][-r random samples][-N Newick formatted tree][-R cut-off R]\n");
				exit(1);
				break;
			case '?':
				fprintf(stderr, "Error: Unknown option passed in: %c\n", optopt);
				fprintf(stderr, "usage: [-F alignmentfolder] options([-S Structure folder]\n[--intra/--inter analysis type][-a threshold alpha value][-r random samples][-N Newick formatted tree][-R cut-off R]\n");
				exit(1);
				break;
		}
	}

	Pwd(output);
	output += "/";

	if(mystring.length()==0){
		fprintf(stderr, "Error: you haven't entered an alignment folder! Use -F option.\n");
		exit(-1);
	}

	vector<string> files;
	files = Folder_to_vector(mystring.c_str());

	Fasta_vector file;
	file.ref_num=0;
	Fasta_map file1, file2;

	/*set up random seed*/
	gettimeofday(&start, NULL);
	seed = start.tv_sec*start.tv_usec;
	gsl_rng_env_setup();
	T=gsl_rng_default;
	r=gsl_rng_alloc(gsl_rng_mt19937);
	gsl_rng_set(r, seed);
	/*finished initialising seed*/
	string B_aa="ACDEFGHIKLMNPQRSTVWY-";
	vector< vector<double> > D,D_corrected; 

	//Blosum Blos62 = Blosum62();
	//Blosum Blos62 = Blosum45();


	TreeTemplate<Node> *fixed = NULL;
	if(tree_in ==1){
		Newick * newickReader = new Newick(false);
		fixed = newickReader->read(treefile.c_str());
		delete newickReader;
	}

	vector<Node *> nams;


	ofstream interact;
	if(analysis==1){
		interact.open("coev_inter.csv");
		interact << "File1\tFile2\tnum_pairs\ttotal_comp\tCut Off\tthreshold r\taverage r\taverage sig r\ttree1length\ttree2length\tgap threshold\tbootcutoff\tDistance Coef\n";
	}
	vector< vector<int> > coev;
	vector< vector<int> > coev_total;
	double aver_data=0.0;

	vector< vector<double> > aver_back;
	vector<double> temp_back(files.size(), 0.0);
	vector<int> temp_ba(files.size(), 0);
	for(unsigned int i=0;i<files.size();++i){
		aver_back.push_back(temp_back);
		coev.push_back(temp_ba);
		coev_total.push_back(temp_ba);
	}


	ofstream intra_nums;
	if(analysis==0){
		intra_nums.open("coev_intra.txt");
		intra_nums << "Filename" << "\tPairs of Sites\t#of comparisons\tthreshold\tboot-cutoff\tTotal tree length" << endl;
	}


	FILE *stream;
	stream = freopen("/dev/null", "w", stdout);

	for(unsigned int i=0;i<files.size();++i){
		output.clear();
		Pwd(output);
		string temp_file;
		temp_file=mystring;
		temp_file+= "/";
		temp_file += files[i];
		vector<double> Correlations;




		if(analysis ==0){
			file.clear();
			Read_Fasta_vector(temp_file.c_str(), file);	
			file.check_for_pdb(pdb_folder, files[i]);
			TreeTemplate<Node> *tree = NULL;

			if(treefol==true){
				string treefile1 = file.check_for_tree(treefile, files[i]);
				Newick * newickReader = new Newick(false);
				string tem;
				Pwd(tem);
				tem += "/";
				tem += treefile1;
				try{
					tree = newickReader->read(tem.c_str());
				}catch(exception & e){
					cout << "No tree file submitted, will create a distance tree\n";
				}
				if(tree==NULL){
					treefol=false;
					cerr << "Couldn't find tree file, will create BioNJ tree" << endl;
				}
				delete newickReader;
			}
			file.ref_num = 0;	


			file.check_valid_alignment("AMINO");
			file.get_identity_level();

			Blosum Blos1;
			if(file.identity<0.475){
				Blos1 = Blosum45();
				//Blos1 = Blosum62();
			}else if(file.identity>0.475&&file.identity<0.575){
				Blos1 = Blosum50();
				//	Blos1 = Blosum62();
			}else if(file.identity>0.575&&file.identity<0.675){
				Blos1 = Blosum62();
			}else if(file.identity>0.675&&file.identity<0.775){
				Blos1 = Blosum70();
				//	Blos1 = Blosum62();
			}else{
				Blos1 = Blosum80();
				//	Blos1 = Blosum62();
			}

			vector<site> pairs_of_sites;
			output += files[i];
			output += ".out";
			cerr << "Input file: " << temp_file << endl;
			cerr << "Output will be sent to: " << output << endl << endl;
			print_splash(output);
			file.print_to_fasta(output.c_str());
			ofstream OUTPUT(output.c_str(), ios::app);
			int length = file.sequences[0].length();			

			/*get the JTT tree that corresponds to the sequences involved*/


			std::auto_ptr<DistanceMatrix> DS;



			if(treefol ==false){
				cerr << "Creating tree...\n\n";
				tree = create_input_tree(file.names, file.sequences);
			}else if(tree_in==1 && variable==1){

				string tempp = TreeTemplateTools::treeToParenthesis(*fixed, false);	
				tree = Remove_names(file.names, fixed);
				tempp = TreeTemplateTools::treeToParenthesis(*tree, false);
				vector<string> newnames = tree->getLeavesNames();
				nams = Put_distances_on_tree(tree, DS);
			}


			vector<double> totaltemp;
			double threshold;

			cerr << "Performing " << num_randoms << " simulations...\n\n";


			double last =10.2;		//just some large number
			double stanlast =10.2;
			for(int j=0;j<num_randoms;++j){
				vector<string> Tnames, tsequences;
				create_seq(tree, Tnames, tsequences, length);

				vector< vector<double> > Dtemp;
				vector<int> temp_diff;
				Diff_aa_column(tsequences, temp_diff, gapthresh);
				Estimate_D(tsequences, Tnames, Dtemp, Blos1, tree);
				vector<double> tempCorrel;
				int numtempcor;

				/*DO INTRA ON THESE*/
				numtempcor = intra(Dtemp, tempCorrel, (int)file.sequences[0].size(), 1);

				totaltemp.insert(totaltemp.begin(), tempCorrel.begin(), tempCorrel.end());
				double mean = average_vec<double>(totaltemp);
				double stan = SD_vf(totaltemp, mean);

				if(converge==1){
					if(fabs(mean-last) <= 0.00001 && fabs(stan-stanlast)<=0.00001){
						break;
					}else{
						if(j==num_randoms-1)
							j=0;
						last = mean;
						stanlast = stan;
					}

				}
				cerr << "\r" << j/num_randoms << "\% finished";

			}

			cerr << "\r100\% finished" << endl;

			cerr << "\nPerforming Analysis...\n\n";
			sort(totaltemp.begin(), totaltemp.end());
			for(unsigned int l=0;l<totaltemp.size();++l){
				if(totaltemp[l]==1.0){
					totaltemp.erase(totaltemp.begin()+l, totaltemp.end());
				}
			}


			int large = floor((totaltemp.size()*(1-threshval))) +1;
			threshold = totaltemp[large];

			/*get the number of columns with an average number of sites greater than the input value*/

			/*simulate a sequence alignment with the same distances as those given and a length equal to that of the input sequences*/
			/*run the intra coevolution analysis on this (random sample times) to get a cut off value*/

			vector< vector<double> > D;
			vector<int> aa_diff;
			Diff_aa_column(file.sequences, aa_diff, gapthresh);
			Estimate_D(file.sequences, file.names, D, Blos1, tree);

			int count=0;
			for(unsigned int h=0;h<aa_diff.size();++h){
				if(aa_diff[h]>0){
					count++;
				}
			}

			vector<double> Correl;
			int numcor;
			numcor = intra(D, Correl, (int)file.sequences[0].size(), 0);


			int ref=0;
			for(unsigned int p=0;p<file.names.size();++p){
				if(file.names[p]==file.ref_seq){
					ref=p;
					break;
				}else if(file.names[p]==refname){
					ref = p;
					break;
				}	

			}


			unsigned int num_pairs=0;
			num_pairs = print_to_file(Correl, D, threshold, file.sequences, aa_diff, OUTPUT, file, file.names[file.ref_num], files[i], file.ref_num, thresholdR, bootcut, numcor);	
			//intra_nums << files[i] << "\t" << num_pairs << "\t" << ((file.sequences[0].size()-file.getAveGaps())*(file.sequences[0].size()-file.getAveGaps()-1)/2) << "\t" << tree->getTotalLength() << endl;
			intra_nums << files[i] << "\t" << num_pairs << "\t" << numcor << "\t" << threshold << "\t" << bootcut << "\t" << tree->getTotalLength() << endl;
			OUTPUT.close();
			if(tree!=NULL)
				delete tree;
		}else{/*===============THIS IS THE INTER SECTION==============*/


			ofstream OUTPUT;
			file1.clear();


			Read_Fasta_map(temp_file.c_str(), file1);	
			/*loop over all files which are not the same as the i-th file*/
			for(unsigned int j=i+1;j<files.size();++j){
				output.clear();
				Pwd(output);
				file2.clear();
				Fasta_vector vec1, vec2;	
				string temp_file2;
				temp_file2 += mystring; 
				temp_file2 += "/";
				temp_file2 += files[j];
				Read_Fasta_map(temp_file2.c_str(), file2);
				map<string, bool> remnames, remnames2;
				Convert_to_vectors(file1, file2, vec1, vec2, remnames, remnames2);	

				if(vec1.names.size()>=4&&vec2.names.size()>=4){

					vec1.check_for_pdb(pdb_folder, files[i]);
					vec2.check_for_pdb(pdb_folder, files[j]);
					/*	if(vec1.has_pdb)
						vec1.get_amino_neighbours(8.0);
						if(vec2.has_pdb)
						vec2.get_amino_neighbours(8.0);
						*/
					TreeTemplate<Node> *tree1 = NULL;
					TreeTemplate<Node> *tree2 = NULL;


					if(treefol==true){
						string treefile1 = vec1.check_for_tree(treefile, files[i]);
						string treefile2 = vec2.check_for_tree(treefile, files[j]);
						Newick * newickReader = new Newick(false);
						Newick * newickReader2 = new Newick(false);
						string tem, tem2;
						Pwd(tem);
						tem += "/";
						tem += treefile1;
						try{
							tree1 = newickReader->read(tem.c_str());
						}catch(exception & e){
							cout << "No tree file submitted, will create a distance tree\n";
						}

						Pwd(tem2);
						tem2 += "/";
						tem2 += treefile2;
						try{
							tree2 = newickReader2->read(tem2.c_str());
						}catch(exception & e){
							cout << "No tree file submitted, will create a distance tree\n";
						}
						if(tree1==NULL||tree2==NULL){
							treefol=false;
							cerr << "Couldn't find tree file, will create BioNJ tree" << endl;
						}else{
							tree_in=1;
							variable=0;
						}




						map<string, bool>::iterator rit;

						string tempp = TreeTemplateTools::treeToParenthesis(*tree1, false);	
						for(rit=remnames.begin();rit!=remnames.end();++rit){
							//		string tempp = TreeTemplateTools::treeToParenthesis(*tree1, false);	
							TreeTemplateTools::dropLeaf(*tree1, rit->first);
						}

						tempp = TreeTemplateTools::treeToParenthesis(*tree1, false);	
						for(rit=remnames2.begin();rit!=remnames2.end();++rit){
							TreeTemplateTools::dropLeaf(*tree2, rit->first);
						}

						delete newickReader;
						delete newickReader2;
					}


					vec1.check_valid_alignment("AMINO");
					vec2.check_valid_alignment("AMINO");
					//	int sim_length1 = vec1.sequences[0].size() -vec1.getAveGaps();
					//	int sim_length2 = vec2.sequences[0].size() -vec2.getAveGaps();



					//vector<double> rel_dist = distances;
					//	sort(rel_dist.begin(), rel_dist.end());


					vector<site> pairs_of_sites;
					output += files[i];
					output += "_"; output += files[j];output += ".out";
					OUTPUT.open(output.c_str(), ios::app);
					cerr << "Input file1: " << files[i] << endl;
					OUTPUT << "Input file1: " << files[i] << endl;
					OUTPUT << "\n\nInput file2: " << files[j] << endl; 
					cerr << "\nInput file2: " << files[j] << endl; 
					cerr << "\nOutput will be sent to: " << output << endl << endl;
					interact << files[i] << "\t" << files[j] << "\t";


					print_splash(output);
					OUTPUT << "\n\File1: " << files[i] << endl;
					vec1.print_to_fasta(output.c_str());
					OUTPUT << "\n\nFile2: " << files[j] << endl;
					vec2.print_to_fasta(output.c_str());
					int length1 = vec1.sequences[0].length();			
					int length2 = vec2.sequences[0].length();

					OUTPUT << "\n\nLength1: " << length1 << endl;
					OUTPUT << "Length2: " << length2 << endl;


					if(tree_in ==0){
						tree1 = create_input_tree(vec1.names, vec1.sequences);
						tree2 = create_input_tree(vec2.names, vec2.sequences);
						
						// Output the CAPS generated trees to the .out file of each pair
						string temptre1 = TreeTemplateTools::treeToParenthesis(*tree1, true);
						string temptre2 = TreeTemplateTools::treeToParenthesis(*tree2, true);
						OUTPUT << "\n" << endl;
						OUTPUT << "CAPS generated tree 1: " << temptre1 << endl;
						OUTPUT << "CAPS generated tree 2: " << temptre2 << endl;
					}/*else if(tree_in ==1 && variable==1){

					   std::auto_ptr<DistanceMatrix> DS;
					   TreeTemplate<Node> *tree1 = NULL;

					   DS=ScoreDist(vec1.names, vec1.sequences);
					   string tempp = TreeTemplateTools::treeToParenthesis(*fixed, false);	
					   tree1 = Remove_names(vec1.names, fixed);
					   nams = Put_distances_on_tree(tree1, DS);



					   tree2 = Remove_names(vec2.names, fixed);
					   nams = Put_distances_on_tree(tree2, DS);

					   }else{
					   tree1 = fixed;
					   tree2 = fixed;
					   }	
					   */

					vec1.get_identity_level();
					vec2.get_identity_level();
					//cerr << "id=" << vec1.identity << endl;
					//cerr << "id=" << vec2.identity << endl;
					Blosum Blos1;
					if(vec1.identity<0.475){
						Blos1 = Blosum45();
						//Blos1 = Blosum62();
					}else if(vec1.identity>0.475&&vec1.identity<0.575){
						Blos1 = Blosum50();
						//Blos1 = Blosum62();
					}else if(vec1.identity>0.575&&vec1.identity<0.675){
						Blos1 = Blosum62();
					}else if(vec1.identity>0.675&&vec1.identity<0.775){
						Blos1 = Blosum70();
						//Blos1 = Blosum62();
					}else{
						Blos1 = Blosum80();
						//	Blos1 = Blosum62();
					}

					Blosum Blos2;
					if(vec2.identity<0.475){
						Blos2 = Blosum45();
						//Blos2 = Blosum62();
					}else if(vec2.identity>0.475&&vec2.identity<0.575){
						Blos2 = Blosum50();
						//Blos2 = Blosum62();
					}else if(vec2.identity>0.575&&vec2.identity<0.675){
						Blos2 = Blosum62();
					}else if(vec2.identity>0.675&&vec2.identity<0.775){
						Blos2 = Blosum70();
						//Blos2 = Blosum62();
					}else{
						Blos2 = Blosum80();
						//Blos2 = Blosum62();
					}
					vector<double> totaltemp;
					vector<int> numbackground;
					double threshold;
					cerr << "Performing " << num_randoms << " simulations...\n\n";
					double last =10.2;
					double stanlast =10.2;
					//					int twice = 0;

					for(int z=0;z<num_randoms;++z){
						vector<string> Tnames1, tsequences1, Tnames2, tsequences2;
						create_seq(tree1, Tnames1, tsequences1, length1);
						create_seq(tree2, Tnames2, tsequences2, length2);
						vector< vector<double> > Dtemp1, Dtemp2;
						vector<int> temp_diff1, temp_diff2;

						if(z%2==0){
							Diff_aa_column(tsequences1, temp_diff1, gapthresh);
							Estimate_D(tsequences1, Tnames1,Dtemp1, Blos1, tree1);
							Diff_aa_column(tsequences2, temp_diff2, gapthresh);
							Estimate_D(tsequences2, Tnames2, Dtemp2, Blos2, tree1);//this is tree1 for a reason!!
						}else{
							Diff_aa_column(tsequences1, temp_diff1, gapthresh);
							Estimate_D(tsequences1, Tnames1,Dtemp1, Blos1, tree2);
							Diff_aa_column(tsequences2, temp_diff2, gapthresh);
							Estimate_D(tsequences2, Tnames2, Dtemp2, Blos2, tree2);//this is tree1 for a reason!!
						}

						vector<double> tempCorrel;
						inter(Dtemp1, Dtemp2, tempCorrel, 1, temp_diff1, temp_diff2);

						totaltemp.insert(totaltemp.begin(), tempCorrel.begin(), tempCorrel.end());

						double mean = average_vec<double>(totaltemp);
						double stan = SD_vf(totaltemp, mean);

						if(converge==1){
							if(fabs(mean-last) <= 0.00001 && fabs(stan-stanlast)<=0.00001){
								//++twice;
								//	if(twice==2){
								break;
								//	}
							}else{
								//	twice=0;
								if(z==num_randoms-1)
									z=0;
								last = mean;
								stanlast = stan;
							}		
						}
					}

					sort(totaltemp.begin(), totaltemp.end());	
					int value = floor(((totaltemp.size())*(1-(threshval))))+1;

					threshold = totaltemp[value];
					totaltempnew = totaltemp;


					/*=======================================================*/
					/*=======================================================*/
					/*=======================================================*/
					cerr << "Performing Analysis...\n\n";
					vector<double> Correl1, Correl2;
					vector< vector<double> > D, D2;
					vector<int> diff1, diff2;
					Diff_aa_column(vec1.sequences, diff1, gapthresh);
					Diff_aa_column(vec2.sequences, diff2, gapthresh);
					Estimate_D(vec1.sequences, vec1.names, D, Blos1, tree1);
					Estimate_D(vec2.sequences, vec2.names, D2, Blos2, tree1);
					inter(D, D2, Correl1, 0, diff1, diff2);

					Estimate_D(vec1.sequences, vec1.names, D, Blos1, tree2);
					Estimate_D(vec2.sequences, vec2.names, D2, Blos2, tree2);
					inter(D, D2, Correl2, 0, diff1, diff2);

					int ref1=0, ref2=0;
					for(unsigned int p=0;p<vec1.names.size();++p){
						if(vec1.names[p]==vec1.ref_seq){
							ref1=p;
							break;
						}else if(vec1.names[p]==refname){
							ref1 = p;
							ref2=p;
							break;
						}	
					}


					for(unsigned int p=0;p<vec2.names.size();++p){
						if(vec2.names[p]==vec2.ref_seq){
							ref2=p;
							break;
						}	
					}

					vector<double> signif;
					int numcor=0;
					double coef = 0.0;

					//int num_pairs = print_inter(Correl, threshold, OUTPUT, vec1.sequences, vec2.sequences, D, D2, diff1, diff2, vec1, vec2, vec1.names[ref1], vec2.names[ref2], files[i], files[j], ref1, thresholdR, signif, bootcut, numcor, coef);
					int num_pairs = print_inter(Correl1, Correl2, threshold, OUTPUT, vec1.sequences, vec2.sequences, D, D2, diff1, diff2, files[i], files[j], ref1, thresholdR, signif, bootcut, numcor);



					//double P_val=1.0;	
					//					int corl = (vec1.sequences[0].size()-vec1.getAveGaps())* (vec2.sequences[0].size()-vec2.getAveGaps());
					//int corl = numcor;
					//int corl = numcor;

					int corl = vec1.sequences[0].size()*vec2.sequences[0].size();

					/*		if(vec1.has_pdb && vec2.has_pdb){
							corl = vec1.pdb.surface_num*vec2.pdb.surface_num;
							}
							*/	
					//					cerr << "corl=" << corl << "\tvec1sur=" << vec1.pdb.surface_num << "\tvec2sur=" << vec2.pdb.surface_num << "\tpairs=" << num_pairs << endl;

					double av_correl = average_vec<double>(Correl1);				
					av_correl += average_vec<double>(Correl2);
					av_correl /= 2;				

					double av_sig = 0.0;
					if(signif.size()>0)
						av_sig = average_vec<double>(signif);				

					interact <<  num_pairs << "\t" << corl << "\t" << threshval << "\t" << threshold << "\t" << av_correl << "\t" << av_sig << "\t" << tree1->getTotalLength() << "\t" << tree2->getTotalLength()  <<  "\t" << gapthresh << "\t" << bootcut <<  "\t" << coef << endl;
					OUTPUT.close();

					//delete DS;
					delete tree2;
					delete tree1;
					//delete fixed;
					aver_data+=(double)num_pairs/(double)Correl1.size();

					coev[i][j]=num_pairs;
					coev[j][i]=num_pairs;
					coev_total[i][j]= corl;
					coev_total[j][i]= corl;

					//	temp_aver.push_back((double)num_pairs/(double)Correl.size());
					aver_back[i][j] = (double)num_pairs/(double)corl;
					aver_back[j][i] = (double)num_pairs/(double)corl;

				}

			}
		}		


	}


	if(analysis==0){
		intra_nums.close();
	}


	if(analysis==1){

		interact << "\n\nWhat does each gene say about the others           \n";
		interact << "=======================================================\n";
                interact <<  "file[i]\tfile[j]\tINTERACT?\tfile_back[i]\texp\tnott\tback_not\tchi_pass\tchi\n";
		
		vector<double> file_back;
		for(unsigned int i=0;i<files.size();++i){

			double av =0.0;
			for(unsigned int j=0;j<files.size();++j){
				av += aver_back[i][j];
			}

			file_back.push_back(av/((double)files.size()-1));

		}

		int k=0;
		for(unsigned int i=0;i<files.size();++i){

			for(unsigned int j=0;j<files.size();++j){

				if(i!=j){
					int nott = coev_total[i][j]-coev[i][j];
					double chi_pass = gsl_cdf_chisq_Pinv(0.95, 1);
					double back_not = (1-file_back[i])*(double)coev_total[i][j];
					double exp = file_back[i]*(double)coev_total[i][j];	
					double chi = pow(((double)coev[i][j]-exp),2)/exp;
					chi += pow(((double)nott-back_not),2)/back_not;

					interact << files[i] << "\t" << files[j];

					if(chi >=gsl_cdf_chisq_Pinv(0.95, 1) && (double)coev[i][j]>(exp)){
						interact << "\tYES" << "\t" << file_back[i] << "\t" << exp << "\t" << nott << "\t" << back_not << "\t" << chi_pass << "\t" << chi << "\n";
					}else if(chi >=gsl_cdf_chisq_Pinv(0.95, 1) && (double)coev[i][j]<(exp)){
						interact << "\tNO" << "\t" << file_back[i] << "\t" << exp << "\t" << nott << "\t" << back_not << "\t" << chi_pass << "\t" << chi << "\n";
					}else{
						interact << "\tNO" << "\t" << file_back[i] << "\t" << exp << "\t" << nott << "\t" << back_not << "\t" << chi_pass << "\t" << chi << "\n";
					}
					if(i!=j)
						++k;

				}

			}
		}


	}	


	interact.close();

	gsl_rng_free(r);
}



/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  Chi_squared
 *  Description:  Deduces the Probability of co-evolution of 2 proteins
 * =====================================================================================
 */
int Chi_squared (int num_pairs, int num_correlations, double background, double& P_val){

	double chi=0.0;

	//	chi +=	(pow((num_pairs-background),2))/background; 
	//	double nott=(num_correlations - num_pairs)/(double)num_correlations;
	//	double back_not = (num_correlations - background)/(double)num_correlations;
	double nott=(num_correlations - num_pairs);
	double back_not = (num_correlations - background);





	/* This is now the G-test, same theoretical assumptions */
	chi += (num_pairs*(log((double)num_pairs/(double)background)));
	chi += (nott*(log((double)nott/(double)back_not)));
	chi *= 2;


	if(chi>= gsl_cdf_chisq_Pinv(0.999, 1)){
		P_val = 0.001;
	}else if(chi>=gsl_cdf_chisq_Pinv(0.99, 1)){
		P_val = 0.01;
	}else if(chi>=gsl_cdf_chisq_Pinv(0.95, 1)){
		P_val = 0.05;
	}


	if(chi >=gsl_cdf_chisq_Pinv(0.95, 1) && num_pairs>(background)){
		return 1;
	}else if(chi >=gsl_cdf_chisq_Pinv(0.95, 1) && num_pairs<(background)){
		return 2;
	}else{
		return 0;
	}


}		/* -----  end of function Chi_squared  ----- */

/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  find_alpha
 *  Description:  Find the index of an element in a vector totaltemp
 *    Help from:  https://www.geeksforgeeks.org/how-to-find-index-of-a-given-element-in-a-vector-in-cpp/
 *                https://stackoverflow.com/questions/8647635/elegant-way-to-find-closest-value-in-a-vector-from-above
 *       Author:  Petar Petrov, University of Turku (Finland); pebope@utu.fi
 * =====================================================================================
 */
double getIndex(std::vector<double> const& v, double K) 
{ 
    auto const it = std::lower_bound(v.begin(), v.end(), fabs(K));
    //auto it = std::upper_bound(v.begin(), v.end(), fabs(K));
  
    if (it != v.end()) { 
        int index = distance(v.begin(), it);
	alphathresh = (((int)1+(double)v.size()-(int)index)/(double)v.size());
	return alphathresh;
        //cerr << index << "\t" << alphathresh << endl; 
    } 
    else { 
        cerr << "ELEMENT NOT FOUND!" << endl; 
    } 
} 




/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  print_inter
 *  Description:  Print the analysis results for inter molecular coevolution
 * =====================================================================================
 */
//int print_inter(vector<double>& Correl, double threshold, ofstream& output, vector<string>& sequences, vector<string>& sequences2, vector< vector<double> >& D, vector< vector<double> >& D2, vector<int>& diff1, vector<int>& diff2, Fasta_vector& file1, Fasta_vector& file2, string& ref_spec1, string& ref_spec2, string& files1, string& files2, int ref, double& thresholdR, vector<double>& signif, double bootcut, int& numcor, double& coef){
int print_inter(vector<double>& Correl1, vector<double>& Correl2, double threshold, ofstream& output, vector<string>& sequences, vector<string>& sequences2, vector< vector<double> >& D, vector< vector<double> >& D2, vector<int>& diff1, vector<int>& diff2, string& files1, string& files2, int ref, double& thresholdR, vector<double>& signif, double bootcut, int& numcor){


	int gaps1=0;
	int cor=0, pairs=0;

	output << endl << endl;

	output << "Coevolving Pairs of amino acid sites\n";
	output << "================================================================================================================================\n";
	output << "Col1(real)\tCol2(real)\tDmean1\t\tDmean2\t\tCorrelation\tBootstrap value\tP-value1\tP-value2\tMean P-value\tCorrelation1\tCorrelation2\n\n";
	output << "================================================================================================================================\n";

	//double mean = average_vec<double>(Correl);
	//double SD = SD_vf(Correl, mean);

	vector< vector<int> > pair;
	map< int, int> gap_map1, gap_map2;

	/*========================================================*/
	/*========================================================*/


	for(unsigned int i=0;i<sequences[0].size();++i){
		if(sequences[ref][i] == '-'){
			++gaps1;
		}
		int gaps2=0;
		if(diff1[i]>0){
			double averDi = average_vec<double>(D[i]);
			for(unsigned int j=0;j<sequences2[0].size();++j){
				if(sequences2[ref][j]=='-'){
					++gaps2;	
				}
				if(diff2[j]>0){
					double averDj = average_vec<double>(D2[j]);
					if((fabs(Correl1[cor])>=threshold && fabs(Correl1[cor])>=thresholdR) && ( fabs(Correl2[cor])>=threshold && fabs(Correl2[cor])>=thresholdR  ) ){

						//	double neighval = Neigh(i-gaps1+1, j+1-gaps2, file1, file2);

						map<int, int>::iterator it1;
						map<int, int>::iterator it2;

						//int re1=0, re2=0;
						/*	if(file1.has_pdb){	
							it1 = file1.pdb.chains[file1.ref_chan_num].res_map.find(i-gaps1+1);
							if(it1!=file1.pdb.chains[file1.ref_chan_num].res_map.end()){	
							int res1 = file1.pdb.chains[file1.ref_chan_num].res_map[i-gaps1+1];
							re1 = file1.pdb.chains[file1.ref_chan_num].aminos[res1].neighbours;
							}else{
							re1=10;
							}
							}
							if(file2.has_pdb){	
							it2 = file2.pdb.chains[file2.ref_chan_num].res_map.find(j-gaps2+1);
							if(it1!=file1.pdb.chains[file1.ref_chan_num].res_map.end()){	
							int res2 = file2.pdb.chains[file2.ref_chan_num].res_map[j-gaps2+1];
							re2 = file2.pdb.chains[file2.ref_chan_num].aminos[res2].neighbours;
							}else{
							re2=10;
							}

							}
							if(re1<=4 && re2 <=4){
							*/
						double bootval = 0.0;
						bootval = Boot(D, D2, i, j, 10000, threshold);

						//	}

						double Alpha1 = getIndex(totaltempnew, Correl1[cor]);
						double Alpha2 = getIndex(totaltempnew, Correl2[cor]);
						//if(bootval>=bootcut && re1<=8 && re2<=8 ){
						if(bootval>=bootcut){
							output << i+1 << "(" << i-gaps1+1 << ")\t\t" << j+1 << "(" << (j+1)-gaps2 << ")\t\t" << averDi << "\t\t" << averDj << "\t" << (Correl1[cor]+Correl2[cor])/2 << "\t" << bootval << "\t" << Alpha1 << "\t" << Alpha2 << "\t" << (Alpha1+Alpha2)/2 << "\t" << Correl1[cor] << "\t" << Correl2[cor] << endl; 
							signif.push_back(((Correl1[cor]+Correl2[cor])/2));
							++pairs;
							vector<int> tem;
							tem.push_back(i+1);
							tem.push_back(j+1);

							gap_map1[i+1]=i-gaps1+1;
							gap_map2[j+1]=j+1-gaps2;

							pair.push_back(tem);
						}
					}
					++numcor;

					}
					//cerr << "i=" << i << "\tj=" << j << "cor=" << cor << endl;
					++cor;
				}
			}
		}
		/*		if(pair.size()>0 && file1.has_pdb==true && file2.has_pdb==true){
				coef = get_dist_coef(pair, file1, file2, gap_map1, gap_map2);
				}
				*/
		//if(pair.size()>0)
		//	non_over_inter(pair, gap_map1, gap_map2, output, ref_spec1, files1, files2);
		if(pair.size()>0)
			group_inter(pair, gap_map1, gap_map2, output, files1, files2);

		return pairs;
	}		/* -----  end of function print_inter  ----- */


	int Neigh(int i, int j, Fasta_vector& fil1, Fasta_vector& fil2){

		if(fil1.has_pdb && fil2.has_pdb){
			int deg1 = fil1.Get_neighbours(12.0, i);	
			int deg2 = fil2.Get_neighbours(12.0, j);	
			if(deg1>=10 || deg2>=10){
				return 1;
			}else{
				return 0;
			}
		}else{
			return 0;
		}



	}


	double get_dist_coef(vector< vector<int> >& pairs, Fasta_vector& file1, Fasta_vector& file2, map< int, int>& gap_map1, map< int, int>& gap_map2){

		double score=0.0;

		vector<int> left, right;

		for(unsigned int i=0;i<pairs.size();++i){
			left.push_back(pairs[i][0]);
			right.push_back(pairs[i][1]);
		}


		int totleft=0;
		for(unsigned int i=0;i<left.size();++i){
			for(unsigned int j=i+1;j<left.size();++j){
				if(get_distance(gap_map1[left[i]], gap_map1[left[j]], file1)<8.0){
					++totleft;
				}
			}
		}



		int totright = 0;
		for(unsigned int i=0;i<right.size();++i){
			for(unsigned int j=i+1;j<right.size();++j){
				if(get_distance(gap_map2[right[i]], gap_map2[right[j]], file2)<8.0){
					++totright;			
				}
			}
		}

		score = (((double)totleft + (double)totright)/((double)pairs.size()*(double)pairs.size()));

		return score;
	}


	double Boot(vector< vector<double> >& D, vector< vector<double> >& D2, int i, int j, int num_boots, double thresh){


		int total=0;

		int siz=D[0].size();

		for(int k=0;k<num_boots;++k){


			vector<double> D_column(siz, 0.0), D_column2(siz, 0.0);
			for(int n=0;n<siz;++n){
				int ran = gsl_rng_uniform_int(r, siz-1);

				D_column[n] = D[i][ran];
				D_column2[n] = D2[j][ran];
			}

			double mean_a = average_vec<double>(D_column);
			double mean_b = average_vec<double>(D_column2);

			if( Correlation(D_column, D_column2, mean_a, mean_b) >=thresh){
				++total;
			}
		}

		//total/=num_boots;
		double boo = (double)total/(double)num_boots;

		return boo;

	}


	void Rem_inner(TreeTemplate<Node>* tree2){
		//bool go = true;
		//			while(go){
		vector<int> nos = tree2->getLeavesId();

		for(unsigned int i=0;i<nos.size();++i){
			if(!tree2->hasNodeName(nos[i])){
				int dar = tree2->getFatherId(nos[i]);
				tree2->getNode(dar)->removeSon(nos[i]);
				i=0;
				nos.clear();nos = tree2->getLeavesId();
			}
		}


		//			}

	}

	void non_over_inter (vector< vector<int> >& pairs, map<int, int >& gap_map1, map<int, int>& gap_map2, ofstream& output, string& files1, string& files2){



		output << "\n\n";
		output << "Groups of coevolving amino acids" << endl;
		output << "================================" << endl << endl;


		map<int, vector<int> > groups_left, groups_right;

		for(unsigned int i=0;i<pairs.size();++i){

			map<int, vector<int> >::iterator git;
			git=groups_left.find(pairs[i][0]);
			if(git==groups_left.end()){
				vector<int> tem;
				tem.push_back(pairs[i][1]);
				groups_left[pairs[i][0]]=tem;
			}else{
				git->second.push_back(pairs[i][1]);
			}

			git=groups_right.find(pairs[i][1]);

			if(git==groups_right.end()){
				vector<int> tem;
				tem.push_back(pairs[i][0]);
				groups_left[pairs[i][1]]=tem;
			}else{
				git->second.push_back(pairs[i][0]);
			}


		}



		output << files1 << " vs " << files2 << "\n===========================================" << endl << endl;

		ofstream fil1;
		string tem; 
		tem.clear();
		tem += files1;
		tem += "_";
		tem += files2;
		tem += ".csv";

		//fil1.open(tem.c_str());

		//fil1 << "sequenceName,msaColumnNumber,residueNumber,pdbNumber";

		map<int, vector<int> >::iterator mit, mit2;
		int p=1;
		//for(mit=groups_left.begin();mit!=groups_left.end();++mit){
		//fil1 << ",Group" << p;
		//++p;	
		//}
		//fil1 << endl;



		int k=1;
		for(mit=groups_left.begin();mit!=groups_left.end();++mit){

			output << "Group " << k << ": " << mit->first << "(" << gap_map1[mit->first] << "):\t";

			for(unsigned int i=0;i<mit->second.size();++i){
				output << mit->second[i] << "(" << gap_map2[mit->second[i]] << ")\t"; 
				//fil1 << ref_spec << "," << mit->second[i] << ",-1,-1";
				p=1;
				for(mit2=groups_left.begin();mit2!=groups_left.end();++mit2){

					if(p==k){
						//fil1 << ",1";
					}else{
						//fil1 << ",-1";
					}	

				}

				//		fil1 << endl;


			}
			++k;
			output << endl;
		}

		//fil1.close();


		output << "File 2 vs File 1" << "===========================================" << endl << endl;

		tem.clear();
		tem += files2;
		tem += "_";
		tem += files1;
		tem += ".csv";


		k=1;
		for(mit=groups_right.begin();mit!=groups_right.end();++mit){


			output << "Group " << k << ": " << mit->first << "(" << gap_map2[mit->first] << "):\t";

			for(unsigned int i=0;i<mit->second.size();++i){
				output << mit->second[i] << "(" << gap_map1[mit->second[i]] << ")\t"; 

				//fil1 << ref_spec << "," << mit->second[i] << ",-1,-1";
				p=1;
				for(mit2=groups_left.begin();mit2!=groups_left.end();++mit2){

					if(p==k){
						//fil1 << ",1";
					}else{
						//fil1 << ",-1";
					}	

				}

				//fil1 << endl;



			}
			output << endl;
			++k;
		}
		output << endl;




	}


	/* 
	 * ===  FUNCTION  ======================================================================
	 *         Name:  group_inter
	 *  Description:  
	 * =====================================================================================
	 */
	void group_inter (vector< vector<int> >& pairs, map<int, int >& gap_map1, map<int, int>& gap_map2, ofstream& output, string& files1, string& files2){


		map<int, vector<int> > groups;
		map<int, vector<int> > groups2;
		for(unsigned int i=0;i<pairs.size();++i){

			map<int, vector<int> >::iterator it;
			it = groups.find(pairs[i][0]);
			if(it==groups.end()){
				vector<int> tem;
				tem.push_back(pairs[i][1]);
				groups[pairs[i][0]]=tem;
			}else{
				it->second.push_back(pairs[i][1]);
			}
		}


		for(unsigned int i=0;i<pairs.size();++i){

			map<int, vector<int> >::iterator it;
			it = groups2.find(pairs[i][1]);
			if(it==groups2.end()){
				vector<int> tem;
				tem.push_back(pairs[i][0]);
				groups2[pairs[i][1]]=tem;
			}else{
				it->second.push_back(pairs[i][0]);
			}
		}

		output << endl << "Overlapping groups of coevolving residues" << endl << "====================================================================" << endl;
		output << "Group #: Coevolving Sites" << endl << endl;

		map<int, vector<int> >::iterator git, git2;


		map< vector<int>, vector<int> > full_map;
		map< vector<int>, vector<int> >::iterator mit;;




		vector< vector<int> > left, right;

		for(git=groups.begin();git!=groups.end();++git){
			vector<int> temp;
			temp.push_back(git->first);
			left.push_back(temp);
			right.push_back(git->second);
		}


		for(git=groups2.begin();git!=groups2.end();++git){
			for(unsigned int i=0;i<left.size();++i){
				vector<int>::iterator vit, vit2;

				vit = find(right[i].begin(), right[i].end(), git->first);
				if(vit!=right[i].end()){
					for(unsigned int j=0;j<git->second.size();++j){
						vit2 = find(left[i].begin(), left[i].end(), git->second[j]);
						if(vit2==left[i].end()){
							vector<int> v(10000);
							vector<int>::iterator it;
							it = set_union (left[i].begin(), left[i].end(), git->second.begin(), git->second.end(), v.begin());
							v.erase(v.begin()+int(it -v.begin()), v.end());

							//	left[i].push_back(git->second[j]);
							left[i].clear();
							left[i]=v;
						}
					}
				}


			}
		}


		for(git=groups.begin();git!=groups.end();++git){
			for(unsigned int i=0;i<left.size();++i){
				vector<int>::iterator vit, vit2;

				vit = find(left[i].begin(), left[i].end(), git->first);
				if(vit!=left[i].end()){
					for(unsigned int j=0;j<git->second.size();++j){
						vit2 = find(right[i].begin(), right[i].end(), git->second[j]);
						if(vit2==right[i].end()){
							vector<int> v(10000);
							vector<int>::iterator it;
							it = set_union (right[i].begin(), right[i].end(), git->second.begin(), git->second.end(), v.begin());

							v.erase(v.begin()+int(it -v.begin()), v.end());
							right[i].clear();
							right[i]=v;

						}

					}
				}


			}
		}


		int k=1;
		map< vector<int>, vector<int> > mymap;
		map< vector<int>, vector<int> >::iterator myit, myit2;


		for(unsigned int i=0;i<right.size();++i){
			sort(left[i].begin(), left[i].end());
			sort(right[i].begin(), right[i].end());
			mymap[left[i]]=right[i];
		}

		int brekkie=1;

		for(myit=mymap.begin();myit!=mymap.end();++myit){
			if(brekkie==1){
				myit=mymap.begin();
				brekkie=0;
			}


			myit2=myit;
			++myit2;
			for(;myit2!=mymap.end();++myit2){
				vector<int> v(10000), v2(10000);
				vector<int>::iterator it, it2;
				it=set_intersection (myit->first.begin(), myit->first.end(), myit2->first.begin(), myit2->first.end(), v.begin());
				it2=set_intersection (myit->second.begin(), myit->second.end(), myit2->second.begin(), myit2->second.end(), v2.begin());

				if(int(it -v.begin())>0 || int(it2 - v2.begin())>0){
					vector<int> templ(20000), tempr(20000);
					v.erase(v.begin()+int(it -v.begin()), v.end());
					v2.erase(v2.begin()+int(it -v2.begin()), v2.end());

					it=set_union(myit->first.begin(), myit->first.end(), myit2->first.begin(), myit2->first.end(), templ.begin());
					it2=set_union(myit->second.begin(), myit->second.end(), myit2->second.begin(), myit2->second.end(), tempr.begin());
					mymap.erase(myit);
					templ.erase(templ.begin()+int(it -templ.begin()), templ.end());
					tempr.erase(tempr.begin()+int(it2 -tempr.begin()), tempr.end());
					mymap.erase(myit2);
					mymap.insert (mymap.begin(), pair<vector<int> ,vector<int> >(templ,tempr));
					myit=mymap.begin();
					brekkie=1;
					break;	
				}
			}

		}



		//	ofstream gro, gro2;
		string tem;
		tem.clear();
		tem += files1;
		tem += "_";
		tem += files2;
		tem += "-overlap.csv";
		//			cerr << "Residue annotations sent to: " << tem << endl << " and : ";
		//	gro.open(tem.c_str());
		tem.clear();
		tem += files2;
		tem += "_";
		tem += files1;
		tem += "-overlap.csv";
		cerr << tem << endl;
		//gro2.open(tem.c_str());



		//gro << "sequenceName,msaColumnNumber,residueNumber,pdbNumber,score,Group\n";
		//gro2 << "sequenceName,msaColumnNumber,residueNumber,pdbNumber,score,Group\n";

		for(myit=mymap.begin();myit!=mymap.end();++myit){

			output << "Group " << k << ": [";
			for(unsigned int i=0;i<myit->first.size();++i){
				output << myit->first[i] << "(" << gap_map1[myit->first[i]] << ")\t";
			}


			output << " | " ;
			for(unsigned int i=0;i<myit->second.size();++i){
				output << myit->second[i] << "(" << gap_map2[myit->second[i]] << ")\t";
			}

			output << " ] \n";
			/*if(file1.has_pdb==true && file2.has_pdb==true){
			  output << "[ ";
			  vector<int> temp = myit->first;
			  output << group_dist_inter(temp, 0, file1, gap_map1, ref_spec, gro, k);
			  output << " | ";

			  output << group_dist_inter(myit->second, 0, file2, gap_map2, ref_spec2, gro2, k);
			  output << " ]\n";

			  }
			  */
			output << "==========================================================\n";

			++k;
		}

		//gro.close();
		//gro2.close();



		/*	int k=1;
			for(git=groups.begin();git!=groups.end();++git){
			output << "Group " << k << ": " << git->first << " [ ";
			vector<int>::iterator vit;
			for(vit=git->second.begin();vit!=git->second.end();++vit){
			output << *vit << "\t";
			}
			output << " ]" << endl;
			++k;

			}	
			*/






	}		/* -----  end of function group_inter  ----- */




	/* 
	 * ===  FUNCTION  ======================================================================
	 *         Name:  Network::add_connection
	 *  Description:  
	 * =====================================================================================
	 */
	void Network::add_connection (int& i, int& j){


		/*find i, if not there add it*/
		map<int, int>::iterator mit1;
		mit1 = node_map.find(i);

		if(mit1==node_map.end()){



			Nod temp;
			vector<int> tem1, tem2;
			tem1.push_back(j);
			temp.node_connect[i]=tem1;
			nodes.push_back(temp);
			node_map[i]=node_map.size()-1;;

		}else{

			vector<int>::iterator vit;

			int tes = mit1->second;
			//vector<int> tem_vec = nodes[tes].node_connect->second;
			//ves = nodes[tes].node_connect->second();	
			vit = find(nodes[tes].node_connect[tes].begin(), nodes[tes].node_connect[tes].end(), j);
			//map<int, int>::iterator mit;
			//mit = node_map.find(j);

			if(vit!=nodes[tes].node_connect[tes].end()){
				nodes[tes].node_connect[tes].push_back(j);
			}

		}

		mit1 = node_map.find(j);
		if(mit1==node_map.end()){
			Nod temp;
			vector<int> tem1, tem2;
			tem2.push_back(i);
			temp.node_connect[i]=tem1;
			nodes.push_back(temp);
			node_map[j]=node_map.size()-1;

		}else{
			vector<int>::iterator vit;
			int tes = mit1->second;
			vit = find(nodes[mit1->second].node_connect[tes].begin(), nodes[mit1->second].node_connect[tes].end(), i);
			if(vit!=nodes[mit1->second].node_connect[tes].end()){
				nodes[mit1->second].node_connect[tes].push_back(i);
			}
		}

	}		/* -----  end of function Network::add_connection  ----- */



	/* 
	 * ===  FUNCTION  ======================================================================
	 *         Name:  Network:print_network
	 *  Description:  
	 * =====================================================================================
	 */
	void Network::print_network (ofstream& file){

		for(unsigned int i=0;i<nodes.size();++i){
			file << "Group" << i << " :";

			map<int, vector<int> >::iterator it;

			for(it=nodes[i].node_connect.begin();it!=nodes[i].node_connect.end();++it){
				file << it->first << "\t";
				for(unsigned int j=0;j<it->second.size();++j){
					file << it->second[j] << "\t";
				}
				file << endl;

			}

		}

	}		/* -----  end of function Network:print_network  ----- */



	/* 
	 * ===  FUNCTION  ======================================================================
	 *         Name:  test_num
	 *  Description:  test how many values are greater than the threshold 
	 * =====================================================================================
	 */
	int test_num (vector<double>& Correl, double thresh){

		int number=0;
		for(int i=Correl.size()-1;i>=0;--i){
			if(fabs(Correl[i])>=fabs(thresh)){
				++number;
			}
		}

		return number;
	}		/* -----  end of function test_num  ----- */




	/* 
	 * ===  FUNCTION  ======================================================================
	 *         Name:  inter
	 *  Description:  Inter Molecular Co-evolution analysis
	 * =====================================================================================
	 */
	void inter(vector< vector<double> >& D, vector< vector<double> >& D2, vector<double>& Correl, int simulate, vector<int>& diff1, vector<int>& diff2){


		for(unsigned int i=0;i<D.size();++i){
			vector<double> D_column(D[0].size(), 0.0), D_column2(D[0].size(), 0.0);
			/*if(correction == 0)
			  {*/
			for(unsigned int m=0;m<D[0].size();++m){
				D_column[m] = D[i][m];
			}
			/*}else
			  {
			  for(unsigned int m=0;m<D[0].size();++m){
			  D_column[m] = D_correct[i][m];
			  }
			  }*/
			double mean_a = average_vec<double>(D_column);
			for(unsigned int j=0;j<D2.size();++j){
				if(diff2[j]==1 && diff1[i]==1){
					/*if(correction == 0)
					  {*/
					for(unsigned int n=0;n<D[0].size();++n){
						D_column2[n] = D2[j][n];
					}
					/*}else
					  {
					  for(unsigned int n=0;n<D[0].size();++n){
					  D_column2[n] = D_correct2[j][n];
					  }
					  }*/

					double mean_b = average_vec<double>(D_column2);
					if(simulate==1){
						Correl.push_back(fabs(Correlation(D_column, D_column2, mean_a, mean_b)));
						//							Correl.push_back(fabs(gsl_stats_correlation(&D_column[0], 1, &D_column2[0], 1,D[0].size()))); 
					}else{
						//gsl_vector_const_view gsl_x = gsl_vector_const_view_array( &D_column[0], D_column.size());
						//gsl_vector_const_view gsl_y = gsl_vector_const_view_array( &D_column2[0], D_column2.size());
						Correl.push_back(Correlation(D_column, D_column2, mean_a, mean_b));
						//							Correl.push_back(gsl_stats_correlation(&D_column[0], 1, &D_column2[0], 1,D[0].size())); 
					}
				}else{
					Correl.push_back(0.0);
				}
			}/*end while j*/
		}/*end while i*/
	}		/* -----  end of function inter  ----- */



	/* 
	 * ===  FUNCTION  ======================================================================
	 *         Name:  check_gaps
	 *  Description:  
	 * =====================================================================================
	 */
	bool check_gaps (vector<string>& seq, vector<string>& seq2, int i)
	{

		bool test=true;
		int count=0;
		for(unsigned int k=0;k<seq.size();++k){
			if(seq[k][i]=='-'){
				count++;
			}
		}

		if((double)count<(0.8*seq.size())){
			test=false;
		}

		if(count>=10){
			test=true;
		}

		count =0;
		for(unsigned int k=0;k<seq2.size();++k){
			if(seq2[k][i]=='-'){
				count++;
			}
		}

		if((double)count<(0.8*seq.size())){
			test=false;
		}

		if(count>=10){
			test=true;
		}

		return test;;
	}		/* -----  end of function check_gaps  ----- */




	/* 
	 * ===  FUNCTION  ======================================================================
	 *         Name:  Convert_to_vectors
	 *  Description:  For inter-coev analysis, converts maps to vectors of common species 
	 * =====================================================================================
	 */
	void Convert_to_vectors(Fasta_map& file1, Fasta_map& file2, Fasta_vector& vec1, Fasta_vector& vec2, map<string, bool>& remnames, map<string, bool>& remnames2){

		map<string, string>::iterator it, fin;

		vec1.filename = file1.filename;
		vec2.filename = file2.filename;
		vec1.has_pdb=false;
		vec2.has_pdb=false;
		for(it=file1.sequences.begin();it!=file1.sequences.end();++it){

			fin = file2.sequences.find(it->first);
			if(fin!=file2.sequences.end()){
				vec1.names.push_back(it->first);
				vec1.sequences.push_back(it->second);
				vec2.names.push_back(fin->first);
				vec2.sequences.push_back(fin->second);
			}else{
				remnames[it->first]=true;
			}
		}

		for(it=file2.sequences.begin();it!=file2.sequences.end();++it){


			fin = file1.sequences.find(it->first);
			if(fin==file1.sequences.end()){
				remnames2[it->first]=true;
			}

		}



	}		/* -----  end of function Convert_to_vectors  ----- */




	/* 
	 * ===  FUNCTION  ======================================================================
	 *         Name:  print_to_file
	 *  Description:  looks for the Correlations, tests if they are greater than threshold
	 * =====================================================================================
	 */

	unsigned int print_to_file (vector<double>& Correl, vector< vector<double> >& D, double threshold, vector<string>& sequences, vector<int>& aa_differences, ofstream& output, Fasta_vector& file, string& ref_spec, string& filesi, int ref, double& thresholdR, double bootcut, int& numcor){
		int gaps1=0, cor=0;




		if(file.has_pdb==true){
			output << "\nCoevolution Analysis\n";
			output << "\n================================================================================================\n";
			output <<  "Position AA1\tPosition AA2\tMean D1\t\tMean D2\t\tCorrelation\tBootstrap value\tPairwise Distance\n";
			output <<  "------------\t------------\t-------\t\t-------\t\t-----------\t---------------\t-----------------\n";
		}else{
			output << "\nCoevolution Analysis\n";
			output << "\n=======================================================================\n";
			output <<  "Position AA1\tPosition AA2\tMean D1\t\tMean D2\t\tCorrelation\tBootstrap value\n";
			output <<  "------------\t------------\t-------\t\t-------\t\t-----------\t---------------\n";
		}

		map<int, int> gap_map;
		numcor=0;
		vector< vector<int> > pairs_of_sites;
		for(unsigned int i=0;i<sequences[0].size()-1;++i){
			if(sequences[ref][i]=='-'){
				++gaps1;
			}

			int gaps2=0;	
			for(unsigned int j=i+1;j<sequences[0].size();++j){
				if(sequences[ref][j]=='-'){
					++gaps2;
				}	

				if(fabs(Correl[cor])>=threshold && aa_differences[i]>0 && aa_differences[j]>0 && fabs(Correl[cor])>=thresholdR){
					double bootval = Boot(D, D, i, j, 10000, threshold);
					if(bootval>=bootcut){
						output << i+1 << "(" << i+1-gaps1 << ")\t" << j+1 << "(" << j+1-(gaps1+gaps2) << ")\t" << average_vec<double>(D[i]) << "\t\t" << average_vec<double>(D[j]) << "\t\t" << Correl[cor] << "\t" << bootval;
						vector<int> tem;
						tem.push_back(i+1);
						tem.push_back(j+1);
						//				output << "\t" << assess_hyd(sequences, i, j) << "\t";
						gap_map[i+1]=i+1-gaps1;
						gap_map[j+1]=j+1-(gaps1+gaps2);

						if(file.has_pdb==true){
							int real_pos1=i-gaps1;
							int real_pos2=j-(gaps1+gaps2);
							output << "\t" << get_distance(real_pos1, real_pos2, file) << "\t";
						}

						output << endl;

						pairs_of_sites.push_back(tem);
					}
					++numcor;
				}else if(aa_differences[i]>0 && aa_differences[j]>0){
					++numcor;
				}
				++cor;
			}

		}


		if(pairs_of_sites.size()>0)
			non_over_intra(pairs_of_sites, output, gap_map, filesi);
		if(pairs_of_sites.size()>0)
			group_intra(pairs_of_sites, output, file, gap_map, filesi, ref_spec);
		//	Hydro(pairs_of_sites, sequences, output);
		return pairs_of_sites.size();

	}		/* -----  end of function print_to_file  ----- */



	/* 
	 * ===  FUNCTION  ======================================================================
	 *         Name:  get_distance
	 *  Description:  
	 * =====================================================================================
	 */
	double  get_distance (int pos1, int pos2, Fasta_vector& file)
	{

		double dist;


		//cerr << "ref_num=" << file.ref_chan_num << endl;
		//cerr << "chainsize=" << file.pdb.chains.size() << endl;
		if(file.pdb.chains[file.ref_chan_num].res_map.find(pos1)==file.pdb.chains[file.ref_chan_num].res_map.end() || file.pdb.chains[file.ref_chan_num].res_map.find(pos2)==file.pdb.chains[file.ref_chan_num].res_map.end()){
			return 999.9;
		}else{
			//cerr << "ref_num=" << file.ref_num << endl;
			//cerr << "pos1=" << pos1 << "\tpos2=" << pos2 << endl;
			//cerr << "meanx1=" << file.pdb.chains[file.ref_chan_num].aminos[file.pdb.chains[file.ref_chan_num].res_map[pos1]].mean_pos.x << "\tmeanx2=" << file.pdb.chains[file.ref_chan_num].aminos[file.pdb.chains[file.ref_chan_num].res_map[pos2]].mean_pos.x << endl; 
			//cerr << "res_map[" << pos1 << "]=" << file.pdb.chains[file.ref_chan_num].res_map[pos1] << "\tres_map[" << pos2 << "]=" << file.pdb.chains[file.ref_chan_num].res_map[pos2] << endl;
			//cerr << "meany1=" << file.pdb.chains[file.ref_chan_num].aminos[file.pdb.chains[file.ref_chan_num].res_map[pos1]].mean_pos.y << "\tmeany2=" << file.pdb.chains[file.ref_chan_num].aminos[file.pdb.chains[file.ref_chan_num].res_map[pos2]].mean_pos.y << endl; 
			//cerr << "meanz1=" << file.pdb.chains[file.ref_chan_num].aminos[file.pdb.chains[file.ref_chan_num].res_map[pos1]].mean_pos.z << "\tmeanz2=" << file.pdb.chains[file.ref_chan_num].aminos[file.pdb.chains[file.ref_chan_num].res_map[pos2]].mean_pos.z << endl; 
			dist = sqrt(
					(pow(file.pdb.chains[file.ref_chan_num].aminos[file.pdb.chains[file.ref_chan_num].res_map[pos1]].mean_pos.x 
					     - file.pdb.chains[file.ref_chan_num].aminos[file.pdb.chains[file.ref_chan_num].res_map[pos2]].mean_pos.x, 2)) 
					+(pow(file.pdb.chains[file.ref_chan_num].aminos[file.pdb.chains[file.ref_chan_num].res_map[pos1]].mean_pos.y 
							-file.pdb.chains[file.ref_chan_num].aminos[file.pdb.chains[file.ref_chan_num].res_map[pos2]].mean_pos.y, 2) )
					+(pow(file.pdb.chains[file.ref_chan_num].aminos[file.pdb.chains[file.ref_chan_num].res_map[pos1]].mean_pos.z 
							-file.pdb.chains[file.ref_chan_num].aminos[file.pdb.chains[file.ref_chan_num].res_map[pos2]].mean_pos.z, 2) ) );

		}
		return dist;
	}		/* -----  end of function get_distance  ----- */



	/* 
	 * ===  FUNCTION  ======================================================================
	 *         Name:  assess_hyd
	 *  Description:  
	 * =====================================================================================
	 */
	string assess_hyd(vector<string>& sequences, int col1, int col2){

		hydro myhyd;
		myhyd.set_hyd();
		string str;


		double meanh=0.0;
		/*first get average hydro and MW*/
		for(unsigned int i=0;i<sequences.size();++i){
			meanh+= myhyd.values[sequences[i][col1]];					
			meanh+= myhyd.values[sequences[i][col2]];					
		}

		meanh /=(sequences.size()*2);

		double max = (155*0.05);
		double lower = meanh - max;
		double upper = meanh + max;
		int t1=0;
		for(unsigned int i=0;i<sequences.size();++i){
			double tempd=myhyd.values[sequences[i][col1]];
			tempd+=myhyd.values[sequences[i][col1]];
			tempd/=2;
			if(tempd<lower || tempd>upper){
				t1=1;
			}

		}
		MW mymw;
		mymw.set_mw();

		meanh=0.0;
		for(unsigned int i=0;i<sequences.size();++i){
			meanh+= mymw.values[sequences[i][col1]];					
			meanh+= mymw.values[sequences[i][col2]];					
		}

		meanh /=(sequences.size()*2);

		max = (129*0.05);
		lower = meanh - max;
		upper = meanh + max;
		int t2=0;


		for(unsigned int i=0;i<sequences.size();++i){
			double tempd=mymw.values[sequences[i][col1]];
			tempd+=mymw.values[sequences[i][col1]];
			tempd/=2;
			if(tempd<lower || tempd>upper){
				t2=1;
			}
		}


		if(t1==0 && t2==0){
			str += "HYD & MW";
		}else if(t1==0){
			str += "HYD";
		}else if(t2==0){
			str += "MW";
		}else{
			str += "NO INFO";
		}
		return str;
	}		/* -----  end of function assess_hyd  ----- */





	/* 
	 * ===  FUNCTION  ======================================================================
	 *         Name:  Hydro
	 *  Description:  
	 * =====================================================================================
	 */
	/*void Hydro (vector< vector<int> >& pairs, vector<string>& sequences, ofstream& output)
	  {
	  hydro myhyd;
	  myhyd.set_hyd();


	  for(int i=0;i<pairs.size();++i){
	//		get_average_hyd		




	}




	}*/		/* -----  end of function Hydro  ----- */



	/* 
	 * ===  FUNCTION  ======================================================================
	 *         Name:  non_over_intra
	 *  Description:  
	 * =====================================================================================
	 */
	void non_over_intra (vector< vector<int> >& pairs, ofstream& output, map<int, int>& gap_map, string& filei){

		map<int, vector<int> > map_group;



		int last = pairs[0][0];
		vector<int> tem_vec;
		tem_vec.push_back(pairs[0][1]);
		map_group[pairs[0][0]] = tem_vec;


		/* make a map of each residue to its coevolving sites*/
		for(unsigned int i=1;i<pairs.size();++i){

			if(pairs[i][0]==last){
				map_group[pairs[i][0]].push_back(pairs[i][1]);
			}else{
				last = pairs[i][0];
				tem_vec.clear();
				tem_vec.push_back(pairs[i][1]);
				map_group[last]=tem_vec;
			}

		}



		map<int, vector<int> >::iterator mit;

		ofstream fil1;
		string tem;
		tem.clear();
		tem += filei;
		tem += ".csv";
		//fil1.open(tem.c_str());

		vector< vector<int> > total_group;


		output << endl << "Groups of coevolving Amino Acids" << endl << "==========================================\n";

		int k=1;

		vector< vector<int> > subs;

		for(mit=map_group.begin();mit!=map_group.end();++mit){
			for(unsigned int i=0;i<mit->second.size();++i){
				vector<int> tem;
				tem.push_back(mit->first);
				tem.push_back(mit->second[i]);
				for(unsigned int j=i+1;j<mit->second.size();++j){
					map<int, vector<int> >::iterator git, git2;
					git = map_group.find(mit->second[i]);
					//	git2 = map_group.find(mit->second[j]);

					if(git!=map_group.end()){
						vector<int>::iterator vit;
						//vit = git->second.find(mit->second[j]);
						vit = find(git->second.begin(), git->second.end(), mit->second[j]);
						if(vit!=git->second.end()){
							tem.push_back(mit->second[j]);
						}
					}	

				}

				if(tem.size()>2 && not_subset(subs, tem)) {
					output << "Group " << k << ":\t";
					for(unsigned int l=0;l<tem.size();++l){

						output << tem[l]  << "(" << gap_map[tem[l]] << ")\t"; 
						//					gro << ref_spec << "," << sites[i] << ",-1,-1," << "," << k << endl; 
					}
					total_group.push_back(tem);
					output << endl;
					output << "=========================================================\n";
					++k;
					subs.push_back(tem);
				}



			}
		}
		output << endl;

		//			cerr << "Residue annotations sent to: " << tem << endl;

		//fil1 << "sequenceName,msaColumnNumber,residueNumber,pdbNumber";
		/*for(unsigned int i=0;i<total_group.size();++i){
		  fil1 << ",Group" << i+1;
		  }	
		  fil1 << endl;



		  for(unsigned int i=0;i<total_group.size();++i){
		  for(unsigned int j=0;j<total_group[i].size();++j){
		  fil1 << ref_spec << "," << total_group[i][j] << ",-1,-1";
		  for(unsigned int l=0;l<total_group.size();++l){
		  if(l==i){
		  fil1 << ",1";
		  }else{
		  fil1 << ",-1";
		  }
		  }
		  fil1 << endl;

		  }
		  }



		  fil1.close();
		  */



	}		/* -----  end of function non_over_intra  ----- */



	/* 
	 * ===  FUNCTION  ======================================================================
	 *         Name:  is_subset
	 *  Description:  
	 * =====================================================================================
	 */
	bool not_subset (vector< vector<int> >& subs, vector<int>& set ){



		for(unsigned int i=0;i<subs.size();++i){
			unsigned int hits=0;

			for(unsigned int j=0;j<set.size();++j){
				if(find(subs[i].begin(), subs[i].end(), set[j])!=subs[i].end()){
					++hits;
				}
			}

			if(hits==set.size()){
				return false;
			}


		}



		return true;
	}		/* -----  end of function is_subset  ----- */



	/* 
	 * ===  FUNCTION  ======================================================================
	 *         Name:  group_intra
	 *  Description:  
	 * =====================================================================================
	 */
	void group_intra (vector< vector<int> >& pairs, ofstream& output, Fasta_vector& file, map<int, int>& gap_map, string& filei, string& ref_spec){


		int cur = pairs[0][0];
		vector< vector<int> > groups;
		vector<int> tg;
		tg.push_back(pairs[0][0]);

		for(unsigned int i=0;i<pairs.size();++i){
			if(pairs[i][0]==cur){
				tg.push_back(pairs[i][1]);		
			}else{
				cur = pairs[i][0];
				groups.push_back(tg);
				tg.clear();
				tg.push_back(pairs[i][0]);
				tg.push_back(pairs[i][1]);

			}


		}


		int reset=0;
		for(unsigned int i=0;i<groups.size();++i){


			for(unsigned int j=i+1;j<groups.size();++j){

				if(reset==1){
					i=0;
					j=1;
					reset=0;
				}

				vector<int> mer(groups[i].size()+groups[j].size());
				merge(groups[i].begin(), groups[i].end(),groups[j].begin(),groups[j].end(), mer.begin());

				map<int, int> mymap;
				for(unsigned int k=0;k<mer.size();++k){
					mymap[mer[k]]=1;
				} 
				if(mymap.size()<(groups[i].size()+groups[j].size())){
					groups.erase(groups.begin()+i);
					groups.erase(groups.begin()+j);
					vector<int> tem;
					map<int, int>::iterator mit;

					for(mit=mymap.begin();mit!=mymap.end();++mit){
						tem.push_back(mit->first);
					}
					groups.push_back(tem);
					reset=1;

				}
			}
		}




		output << endl << endl << "Overlapping groups of coevolving amino acids" << endl << "=================================================" << endl;

		ofstream gro;
		string tem;
		tem += filei;
		tem += "-overlap.csv";
		//			cerr << "Over-lapping residue annotations sent to: " << tem << endl;

		//			gro.open(tem.c_str());
		//			gro << "sequenceName,msaColumnNumber,residueNumber,pdbNumber,score,Group\n";

		for(unsigned int i=0;i<groups.size();++i){
			output << "Group " << i+1 << ": ";
			for(unsigned int j=0;j<groups[i].size();++j){
				output << groups[i][j] << "(" << gap_map[groups[i][j]] << ")\t"; 
			}

			if(file.has_pdb==true){
				output << "mean distance: " << group_dist_intra(groups[i], 0, file, gap_map, ref_spec);
			}


			output << endl;
			output << "==================================================" << endl;
		}

		//gro.close();



	}		/* -----  end of function group_intra  ----- */


	/* 
	 * ===  FUNCTION  ======================================================================
	 *         Name:  group_dist
	 *  Description:  
	 * =====================================================================================
	 */
	double group_dist_intra (vector<int>& sites, int site, Fasta_vector& file , map<int, int>& gap_map, string& species)
	{
		double mean=0.0;
		int k=0;

		for(unsigned int i=0;i<sites.size();++i){
			for(unsigned int j=i+1;j<sites.size();++j){	
				if(check_existance(gap_map[sites[i]], gap_map[sites[j]], file)){
					double temp = get_distance(gap_map[sites[i]], gap_map[sites[j]], file);

					if(temp!=999.9){
						mean += temp;
						++k;
					}
				}
			}
		}
		if(site!=0){
			for(unsigned int i=0;i<sites.size();++i){
				double temp = get_distance(gap_map[site], gap_map[sites[i]], file);

				if(temp!=999.9){
					mean += temp;
					++k;
				}
			}
		}
		if(k!=0 && sites.size()>=1){
			mean /=k;
		}else{
			mean =999.9;
		}

		string temse;
		temse += species;
		for(unsigned int i=0;i<sites.size();++i){
			if(check_existance(gap_map[sites[i]], gap_map[sites[i]], file)){

				//gro << temse << "," << sites[i] << ",-1,-1," << mean << "," << group_num+1 << endl; 
			}
		}
		return mean;
	}		/* -----  end of function group_dist  ----- */


	/* 
	 * ===  FUNCTION  ======================================================================
	 *         Name:  group_dist_inter
	 *  Description:  
	 * =====================================================================================
	 */
	double group_dist_inter (vector<int>& sites, int site, Fasta_vector& file , map<int, int>& gap_map, string& species)
	{
		double mean=0.0;
		int k=0;



		for(unsigned int i=0;i<sites.size();++i){
			for(unsigned int j=i+1;j<sites.size();++j){	
				if(check_existance(gap_map[sites[i]], gap_map[sites[j]], file)){
					double temp = get_distance(gap_map[sites[i]], gap_map[sites[j]], file);

					if(temp!=999.9){
						mean += temp;
						++k;
					}
				}
			}
		}

		if(site!=0){
			for(unsigned int i=0;i<sites.size();++i){
				double temp = get_distance(gap_map[site], gap_map[sites[i]], file);

				if(temp!=999.9){
					mean += temp;
					++k;
				}
			}
		}
		if(k!=0 && sites.size()>=1){
			mean /=k;
		}else{
			mean =999.9;
		}

		for(unsigned int i=0;i<sites.size();++i){
			if(check_existance(gap_map[sites[i]], gap_map[sites[i]], file)){

				string temse;
				temse += species;
				//gro << temse << "," << sites[i] << ",-1,-1," << mean << "," << group_num << endl; 
			}
		}
		return mean;
	}		/* -----  end of function group_dist_inter  ----- */


	/* 
	 * ===  FUNCTION  ======================================================================
	 *         Name:  check_existance
	 *  Description:  
	 * =====================================================================================
	 */
	bool check_existance (int site, int site2, Fasta_vector& file)
	{


		if(file.pdb.chains[0].res_map.find(site)!=file.pdb.chains[0].res_map.end() && file.pdb.chains[0].res_map.find(site2)!=file.pdb.chains[0].res_map.end()){

			return true;
		}else{
			return false;
		}


	}		/* -----  end of function check_existance  ----- */

	/* 
	 * ===  FUNCTION  ======================================================================
	 *         Name:  comparison_fabs
	 *  Description:  Absolute value comparison
	 * =====================================================================================
	 */

	bool comparison_fabs(double i,double j){ 
		return (fabs(i)<fabs(j));
	}


	/* -----  end of function comparison_fabs  ----- */
	/* 
	 * ===  FUNCTION  ======================================================================
	 *         Name:  SD_vf
	 *  Description:  Standard deviation of a vector of doubles
	 * =====================================================================================
	 */


	double SD_vf(vector<double>& array, double mean){
		double total=0.0;
		int size = array.size();

		for(int i=0;i<size;++i){

			total += pow(((double)array[i] - mean),2);

		}
		total =total/size;

		return sqrt(total);
	}

	/* -----  end of function SD_vf  ----- */


	/* 
	 * ===  FUNCTION  ======================================================================
	 *         Name:  Diff_aa_column
	 *  Description:  find differences in the columns
	 * =====================================================================================
	 */


	int Diff_aa_column(vector<string>& sequences, vector<int>& array, double gapthresh){;

		int seq_number = sequences.size();
		int length_sequences = sequences[0].size();
		array.assign(length_sequences, 0);

		vector<int> gapsval;
		for(int i=0;i<length_sequences;++i){

			int totgap = 0;
			for(int j=0;j < seq_number;++j){
				if(sequences[j][i]=='-'){
					totgap++;
				}	

			}
			if((1-(totgap/seq_number))>=gapthresh){
				gapsval.push_back(0);
			}else{
				gapsval.push_back(1);
			}

		}


		for(int i=0;i<length_sequences;++i){

			if(gapsval[i]==0){
				map<char, int> diffmap;

				for(int j=0;j < seq_number;++j){
					map<char, int>::iterator it;
					it = diffmap.find((char)sequences[j][i]);
					if(it==diffmap.end()){
						if(sequences[j][i]!='-')
							diffmap[(char)sequences[j][i]]=1;
					}else{
						++diffmap[(char)sequences[j][i]];
					}


				}//end for j


				if(diffmap.size()>3){
					array[i]=1;
				}else if(diffmap.size()>1){
					int total=0;
					map<char, int>::iterator mit;
					for(mit=diffmap.begin();mit!=diffmap.end();++mit){
						if(mit->second>=2){
							++total;
						}
					}
					if(total>=2){
						array[i]=1;
						//	break;

					}
				}
			}

		}//end i
		return 0;
	}


	/* -----  end of function Diff_aa_column  ----- */


	/* 
	 * ===  FUNCTION  ======================================================================
	 *         Name:  D_average
	 *  Description:  Gets the average of each row of D
	 * =====================================================================================
	 */


	int D_average(vector< vector<double> >& D, vector<double>& D_aver){
		int i=0, j=0, totals;

		totals = D[0].size();
		int size_D = D.size();

		D_aver.assign(size_D, 0.0);

		while(i<size_D-1){
			j=0;
			while(j<=totals){
				D_aver[i] += D[i][j];
				++j;

			}
			D_aver[i]/=totals;
			++i;
		}
		return 0;
	}

	/* -----  end of function D_average  ----- */


	/* 
	 * ===  FUNCTION  ======================================================================
	 *         Name:  average_vi
	 *  Description:  Gives average of a vector of ints, returns double
	 * =====================================================================================
	 */

	double average_vi(vector<int>& array){
		double total;
		int size = array.size();

		total=0.0;
		for(int i=0; i<size;++i){

			total += array[i];
		}
		total /= size;
		return total; 
	}
	/* -----  end of function average_vi  ----- */



	/* 
	 * ===  FUNCTION  ======================================================================
	 *         Name:  getAncestralSequences
	 *  Description:  
	 * =====================================================================================
	 */



	int getAncestralSequences(int args, char ** argv, TreeTemplate<Node>* mytree, vector<string>& sequences, vector<string>& names, map<int, string>& mapinner)
	{

		try {

			BppApplication bppancestor(args, argv, "BppAncestor");
			bppancestor.startTimer();
			Alphabet* alphabet = SequenceApplicationTools::getAlphabet(bppancestor.getParams(), "", false);
			VectorSiteContainer *allSites = new VectorSiteContainer(alphabet);

			SequenceContainer *myseq = SequenceContainerTools::createContainerOfSpecifiedSize(alphabet, sequences.size());

			for(unsigned int i=0;i<sequences.size();++i){
				Sequence temseq(names[i], sequences[i], alphabet);
				myseq->addSequence(temseq, true);
				allSites->addSequence(temseq, alphabet);
			}	
			delete myseq;


			VectorSiteContainer* sites = SequenceApplicationTools::getSitesToAnalyse(*allSites, bppancestor.getParams());
			delete allSites;

			// Get the initial tree
			Tree* tree = mytree->clone();
			string treeWIdPath = ApplicationTools::getAFilePath("output.tree_ids.file", bppancestor.getParams(), false, false);
			DRTreeLikelihood *tl;
			string nhOpt = ApplicationTools::getStringParameter("nonhomogeneous", bppancestor.getParams(), "general", "", true, false);
			ApplicationTools::displayResult("Heterogeneous model", nhOpt);

			SubstitutionModel    *model    = 0;
			model = new JTT92(dynamic_cast<ProteicAlphabet *>(alphabet)); 
			SubstitutionModelSet *modelSet = 0;
			DiscreteDistribution *rDist    = 0;
			unsigned int nbStates;

			FullFrequenciesSet* fSet = new FullFrequenciesSet(model->getAlphabet());
			fSet->setNamespace("anc.");
			fSet->setFrequencies(model->getFrequencies());
			modelSet = SubstitutionModelSetTools::createHomogeneousModelSet(dynamic_cast<SubstitutionModel*>(model->clone()), fSet, tree);
			rDist = PhylogeneticsApplicationTools::getRateDistribution(bppancestor.getParams());
			tl = new DRNonHomogeneousTreeLikelihood(*tree, *sites, modelSet, rDist, true);
			nbStates = modelSet->getNumberOfStates();

			tl->initialize();

			delete tree;
			tree = new TreeTemplate<Node>(tl->getTree());
			bool probs = false;

			AncestralStateReconstruction *asr = NULL;
			asr = new MarginalAncestralStateReconstruction(tl);
			probs = ApplicationTools::getBooleanParameter("asr.probabilities", bppancestor.getParams(), false, "", true, false);

			// Write infos to file:
			string outputFile = ApplicationTools::getAFilePath("output.sites.file", bppancestor.getParams(), false, false);
			outputFile = ApplicationTools::getAFilePath("output.nodes.file", bppancestor.getParams(), false, false);


			SequenceContainer *asSites = 0;
			asSites = asr->getAncestralSequences();

			vector<int> inner;
			TreeTemplate<Node> trel(*tree);
			inner = TreeTemplateTools::getInnerNodesId(*trel.getNode(trel.getRootId()));

			for(unsigned int i=0;i<inner.size();++i){
				Sequence  *my = asr->getAncestralSequenceForNode(inner[i]);
				mapinner[inner[i]] = my->toString(); 
				delete my;
			}


			delete asSites;

			delete asr;
			delete alphabet;
			delete sites;
			if(model)    delete model;
			if(modelSet) delete modelSet;
			delete rDist;
			delete tl;
			delete tree;
			bppancestor.done();

		}
		catch(exception & e)
		{
			cout << e.what() << endl;
			return 1;
		}

		return 0;
	}

	/* -----  end of function getAncestralSequences  ----- */



	void Estimate_D_inter(vector<string>& seqes, vector<string>& names, vector< vector<double> >& D, Blosum& Blos62, TreeTemplate<Node>* tree){

		int total_theta;

		char **arg;


		int ar = 9;

		arg = (char**)malloc(ar*sizeof(char *));

		for(int i=0;i<ar;++i){
			arg[i] = (char*)malloc(256*sizeof(char));
		}

		sprintf(arg[0], "%s", "bppancestor");
		sprintf(arg[1], "%s", "input.sequence.file=fakefile");
		sprintf(arg[2], "%s", "alphabet=Protein");
		sprintf(arg[3], "%s", "input.sequence.sites_to_use=all");
		sprintf(arg[4], "%s", "input.sequence.max_gap_allowed=100%");
		sprintf(arg[5], "%s", "input.tree.file=yeast_tree.tre");
		sprintf(arg[6], "%s", "model=JTT92");
		sprintf(arg[7], "%s", "output.sequence.file=test2");
		sprintf(arg[8], "%s", "nonhomogenous=general");


		map<int, string> mapinner;
		getAncestralSequences(ar, arg, tree, seqes, names, mapinner);



		for(int i=0;i<ar;++i){
			free(arg[i]);
		}
		free(arg);



		int num_seqs = seqes.size();
		int length_seqs = seqes[0].size();
		total_theta = num_seqs*(num_seqs-1);
		//	int number_comp = (num_seqs*(num_seqs -1))/2;

		vector<int> nodes = tree->getNodesId();
		int number_comp  = nodes.size()*nodes.size();;

		vector<double> tempD(number_comp, 0.0);
		vector<int> leaves = tree->getLeavesId();
		map<string, int> names_id;

		for(unsigned int i=0;i<leaves.size();++i){
			names_id[tree->getNodeName(leaves[i])] = leaves[i];
		}

		for(int i=0;i<length_seqs;++i){
			D.push_back(tempD);

			//int k=0;
			vector<int> theta;
			vector<int> fathers;
			//int test=1;
			//	cerr << "line=" << __LINE__ << endl;
			for(unsigned int j=0;j<names.size();++j){

				/*int nod = names_id[names[j]]; 
				  int darth = tree->getFatherId(nod);
				  if(darth!=tree->getRootId()){
				  string temp, temp2;
				  temp += mapinner[darth][i];
				  temp2 += seqes[j][i];
				  theta.push_back(Blos62.values[Blos62.index[temp2]][Blos62.index[temp]]);
				  }

*/

				for(unsigned int k=0;k<names.size();++k){
					string temp, temp2;
					temp += seqes[k][i];
					temp2 += seqes[j][i];
					theta.push_back(Blos62.values[Blos62.index[temp2]][Blos62.index[temp]]);
					//theta.push_back(Blos62.values[Blos62.index[seqes[j][i].c_str()]][Blos62.index[seqes[k][i]]])
				}


				map<int, string>::iterator mit;

				for(mit=mapinner.begin();mit!=mapinner.end();++mit){
					string temp, temp2;
					temp += mapinner[mit->first][i];
					temp2 += seqes[j][i];
					theta.push_back(Blos62.values[Blos62.index[temp2]][Blos62.index[temp]]);
					//theta.push_back(Blos62.values[Blos62.index[mapinner[mit->first][i]]][Blos62.index[seqes[k][i]]])
				}

			}

			//	cerr << "line=" << __LINE__ << endl;
			/*				map<int, string>::iterator mit;

							for(mit=mapinner.begin();mit!=mapinner.end();++mit){

							string temp, temp2;
							int da = tree->getFatherId(mit->first);	
							temp += mapinner[da][i];
							temp2  += mapinner[mit->first][i];
							theta.push_back(Blos62.values[Blos62.index[temp2]][Blos62.index[temp]]);	
							}


			//sequences.push_back(



			}*/


			double Theta_av = average_vi(theta);

			for(unsigned int j=0;j<theta.size();++j){
				D[i][j] = ((theta[j] - Theta_av)*(theta[j] - Theta_av));
				//				D[i][j] = theta[j];
			}

	}

}


/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  Estimate_D
 *  Description:  Estimation of D
 * =====================================================================================
 */


void Estimate_D(vector<string>& seqes, vector<string>& names, vector< vector<double> >& D, Blosum& Blos62, TreeTemplate<Node>* tree){

	int total_theta;

	char **arg;


	int ar = 9;

	arg = (char**)malloc(ar*sizeof(char *));

	for(int i=0;i<ar;++i){
		arg[i] = (char*)malloc(256*sizeof(char));
	}

	sprintf(arg[0], "%s", "bppancestor");
	sprintf(arg[1], "%s", "input.sequence.file=fakefile");
	sprintf(arg[2], "%s", "alphabet=Protein");
	sprintf(arg[3], "%s", "input.sequence.sites_to_use=all");
	sprintf(arg[4], "%s", "input.sequence.max_gap_allowed=100%");
	sprintf(arg[5], "%s", "input.tree.file=yeast_tree.tre");
	sprintf(arg[6], "%s", "model=JTT92");
	sprintf(arg[7], "%s", "output.sequence.file=test2");
	sprintf(arg[8], "%s", "nonhomogenous=general");


	map<int, string> mapinner;
	getAncestralSequences(ar, arg, tree, seqes, names, mapinner);



	for(int i=0;i<ar;++i){
		free(arg[i]);
	}
	free(arg);

	//	cerr << "line=" << __LINE__ << endl;


	int num_seqs = seqes.size();
	int length_seqs = seqes[0].size();
	total_theta = num_seqs*(num_seqs-1);
	//	int number_comp = (num_seqs*(num_seqs -1))/2;

	//	cerr << "line=" << __LINE__ << endl;
	vector<int> nodes = tree->getNodesId();
	int number_comp  = nodes.size()-1;

	vector<double> tempD(number_comp, 0.0);
	vector<int> leaves = tree->getLeavesId();
	map<string, int> names_id;

	//	cerr << "line=" << __LINE__ << endl;
	for(unsigned int i=0;i<leaves.size();++i){
		names_id[tree->getNodeName(leaves[i])] = leaves[i];
	}

	//	cerr << "line=" << __LINE__ << endl;
	for(int i=0;i<length_seqs;++i){
		D.push_back(tempD);

		//int k=0;
		vector<int> theta;
		vector<int> fathers;
		//int test=1;
		//	cerr << "line=" << __LINE__ << endl;
		for(unsigned int j=0;j<names.size();++j){

			//cerr << "nam[" << j << "]=" << names[j] << endl;

			/*		map<string, int>::iterator kit;
					for(kit=names_id.begin();kit!=names_id.end();++kit){
					cerr << "kinam=" << kit->first << "\tsec=" << kit->second << endl;
					}
					*/
			//	cerr << "line=" << __LINE__ << endl;
			int nod = names_id[names[j]]; 
			////cerr << "names[" << j << "]=" << names[j] << endl;
			//	cerr << "line=" << __LINE__ << endl;
			int darth = tree->getFatherId(nod);
			//	cerr << "line=" << __LINE__ << endl;
			if(darth!=tree->getRootId()){
				//		fathers.push_back(darth);
				string temp, temp2;
				temp += mapinner[darth][i];
				temp2 += seqes[j][i];
				//			int test = Blos62.values[Blos62.index[temp2]][Blos62.index[temp]];
				theta.push_back(Blos62.values[Blos62.index[temp2]][Blos62.index[temp]]);
				//				theta.push_back(test);
			}
		}

		map<int, string>::iterator mit;

		for(mit=mapinner.begin();mit!=mapinner.end();++mit){

			if(mit->first!=tree->getRootId()){
				string temp, temp2;
				int da = tree->getFatherId(mit->first);	
				temp += mapinner[da][i];
				temp2  += mapinner[mit->first][i];
				theta.push_back(Blos62.values[Blos62.index[temp2]][Blos62.index[temp]]);	
			}

		}


		double Theta_av = average_vi(theta);

		for(unsigned int j=0;j<theta.size();++j){
			D[i][j] = ((theta[j] - Theta_av)*(theta[j] - Theta_av));
			//				D[i][j] = theta[j];
		}

	}

}

/* -----  end of function Estimate_D  ----- */

/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  Correlation
 *  Description:  Pearsons R
 * =====================================================================================
 */

double Correlation(vector<double>& sample_a, vector<double>& sample_b, double& mean_a, double& mean_b){
	double nominator=0.0, X=0.0, Y=0.0, correl=0.0;


	for(unsigned int i=0;i<sample_a.size();++i){
		nominator += (sample_a[i] - mean_a) * (sample_b[i] - mean_b);
		Y += (sample_b[i] - mean_b)* (sample_b[i] - mean_b);
		X += (sample_a[i] - mean_a)*(sample_a[i] - mean_a);
	}

	//		test = sqrt(X*Y);
	//if(test!=0.0){
	if(X!=0.0 && Y!=0.0){
		//correl = (nominator)/test;
		correl = (nominator)/(sqrt(X*Y));
	}
	return correl;
}

/* -----  end of function Correlation  ----- */

/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  intra
 *  Description:  Intra-molecular co-evolution
 * =====================================================================================
 */


double intra(vector< vector<double> >& D, vector<double>& Correl, int length_seq, int simulate){
	unsigned int j=1;
	int m=0, n=0, number_correlations=0, size_D_column, total_number_comp;
	double percent;

	size_D_column = D[0].size();

	total_number_comp =(( (length_seq) - 1) * length_seq)/2;
	vector<double> D_column(size_D_column, 0.0), D_column2(size_D_column, 0.0);



	/* change these loops so that we have the whole i-loop twice within each of if and else, profiling!! */
	for(unsigned int i=0;i<D.size()-1;++i){
		m = 0;

		//	if(information13 == 0){
		while(m <= size_D_column - 1){
			D_column[m] = D[i][m];
			m++;
		}

		/*	}else{

			while(m <= size_D_column - 1){
			D_column[m] = D_correct[i][m];
			m++;
			}
			}
			*/
		double mean_a = average_vec<double>(D_column);
		for(j=i+1;j< D.size();++j){
			n = 0;

			//		if(information13 == 0){
			while(n <= size_D_column - 1){
				D_column2[n] = D[j][n];
				n++;
			}
			/*	}else{
				while(n <= size_D_column - 1){
				D_column2[n] = D_correct[j][n];
				n++;
				}
				}
				*/
			double mean_b = average_vec<double>(D_column2);
			if(simulate==1){
				Correl.push_back(fabs(Correlation(D_column, D_column2, mean_a, mean_b)));
			}else{
				Correl.push_back(Correlation(D_column, D_column2, mean_a, mean_b));
			}

			++number_correlations;
		}
		percent = ((double)i)/ ((double)(length_seq-1)) * 100;
		j = i + 1;
	}

	return number_correlations;
}


/* -----  end of function intra  ----- */




/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  create_input_tree
 *  Description:  makes an input tree which can be scored
 * =====================================================================================
 */

TreeTemplate<Node> *create_input_tree(vector<string>& seq_names, vector< string >& sequences){
	FILE *stream;
	stream = freopen("/dev/null", "w", stdout);

	std::auto_ptr<DistanceMatrix> DS;
	vector<string> names;

	int tot=0;
	names.resize(seq_names.size());
	for(int i=0;i<(signed)seq_names.size();++i){
		names[i] = seq_names[i];
		++tot;
	}

	DS=ScoreDist(names, sequences);
	/*if(DS==NULL){
	  return NULL;
	  }*/

	AgglomerativeDistanceMethod * distMethod = NULL;
	BioNJ * bionj = new BioNJ();
	bionj->outputPositiveLengths(true);
	distMethod = bionj;

	bionj->setDistanceMatrix(*DS);
	bionj->computeTree(true);
	//if(DS!=NULL)
	//	delete DS;
	TreeTemplate<Node> * tree2 = dynamic_cast<TreeTemplate<Node> *>(bionj->getTree());

	delete bionj;
	return tree2;
}


/*void read_hydros(char *filename){
  fstream hyd;

  hyd.open(filename);





  }*/

int print_splash(string filename){


	if(filename=="stdout"){

	/*	fprintf(stdout, "\n\t\t*******************************************************************************************\n\t\t");
		fprintf(stdout, "* CAPS: Co-Evolution Analysis using Protein Sequences                                     *\n\t\t");
		fprintf(stdout, "* Author: Brian E. Caffrey								  *\n\t\t");
		fprintf(stdout, "* Evolutionary Genetics and Bioinformatics Laboratory					  *\n\t\t");
		fprintf(stdout, "* Smurfit Institute of Genetics								  *\n\t\t");
		fprintf(stdout, "* University of Dublin, Trinity College							  *\n\t\t");
		fprintf(stdout, "* Mathematical Model: Caffrey, Hokamp, Williams and Fares, 2012		  			  *\n\t\t");
		fprintf(stdout, "*******************************************************************************************\n");
	*/	
	}else{

		ofstream OUTPUT(filename.c_str());
		OUTPUT << "\n\t\t*******************************************************************************************\n\t\t";
		OUTPUT << "* CAPS: Co-Evolution Analysis using Protein Sequences                                     *\n\t\t";
		OUTPUT << "* Author: Brian E. Caffrey									  *\n\t\t";
		OUTPUT << "* Evolutionary Genetics and Bioinformatics Laboratory					  *\n\t\t";
		OUTPUT << "* Smurfit Institute of Genetics								  *\n\t\t";
		OUTPUT << "* University of Dublin, Trinity College							  *\n\t\t";
		OUTPUT << "* Mathematical Model: Caffrey, Hokamp, Williams and Fares, 2012  					*\n\t\t";
		OUTPUT << "*******************************************************************************************\n";
	}

	return 0;
}


/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  Remove_names
 *  Description:  
 * =====================================================================================
 */
TreeTemplate<Node>* Remove_names (vector<string>& names, TreeTemplate<Node>* tree){

	vector<int> nodes = tree->getNodesId();

	for(unsigned int i=0;i<nodes.size();++i){
		tree->deleteDistanceToFather(nodes[i]);
	}


	vector<string> treenames = tree->getLeavesNames();
	vector<int> indices;

	/*get indices we don't want*/

	for(unsigned int i=0;i<treenames.size();++i){
		int found =0;
		for(unsigned int j=0;j<names.size();++j){

			if(names[j]==treenames[i]){
				found =1;
				break;
			}
		}

		if(found==0){
			indices.push_back(i);
		}

	}


	string treestring = TreeTemplateTools::treeToParenthesis(*tree, false);

	for(unsigned int i=0;i<indices.size();++i){

		size_t found;
		found = treestring.find(treenames[indices[i]]);

		if(found!=string::npos){
			int start = int(found);
			int end = int(found) + treenames[indices[i]].size();

			string left = whats_left(start, treestring);
			string right = whats_right(end, treestring);

			if(left=="comma" && right=="bracket"){
				int pos=0;
				find_left_bracket(pos, start, treestring);
				treestring.erase(start+1, end-start);
				treestring.erase(pos, 1);

			}else if(left=="bracket" && right=="comma"){
				int pos=0;
				find_right_bracket(pos, end, treestring);
				treestring.erase(pos, 1);
				treestring.erase(start, end-start+1);

			}else if(left=="bracket" && right=="bracket"){
				int pos=0;
				find_comma(pos, start, end, treestring);
				if(pos<start){
					treestring.erase(start, end-start);
					treestring.erase(pos, 1);
				}else{
					treestring.erase(pos, 1);
					treestring.erase(start, end-start);
				}
			}

		}

	}

	TreeTemplate<Node> *newtree;
	stringstream strstr(treestring);

	Newick * newickReader = new Newick(false);
	newtree = newickReader->read(strstr);
	//free(strstr);
	delete newickReader;

	return newtree;
}		/* -----  end of function Remove_names  ----- */


/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  find_comma
 *  Description:  
 * =====================================================================================
 */
void find_comma (int& pos, int start, int end, string tree){

	int left=0, right=0;

	for(int i=start;i>=0;--i){
		if(tree[i]==','){
			left = start - i;
			break;
		}
	}

	for(unsigned int i=end;i<tree.size();++i){
		if(tree[i]==','){
			right =  i - end;
			break;
		}
	}

	if(left<right){
		pos = start - left;
	}else{
		pos = end + right;
	}

}		/* -----  end of function find_comma  ----- */

/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  find_right_bracket
 *  Description:  
 * =====================================================================================
 */
void find_right_bracket (int& pos, int end, string tree){

	int right=0, left=0;

	for(unsigned int i=end;i<tree.size();++i){

		if(tree[i]==')'){
			++right;
		}
		if(tree[i]=='('){
			++left;
		}
		if(right>left){
			pos = i;
			return;
		}
	}


}		/* -----  end of function find_left_bracket  ----- */


/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  find_left_bracket
 *  Description:  
 * =====================================================================================
 */
void find_left_bracket (int& pos, int start, string tree){

	int right=0, left=0;

	for(int i=start;i>=0;--i){

		if(tree[i]==')'){
			++right;
		}
		if(tree[i]=='('){
			++left;
		}
		if(left>right){
			pos = i;
			return;
		}
	}


}		/* -----  end of function find_left_bracket  ----- */


/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  whats_right
 *  Description:  
 * =====================================================================================
 */
string whats_right (int& end, string treestring){

	for(unsigned int i=end;i<treestring.size();++i){

		if(treestring[i]==','){
			end = i;
			return (string)"comma";
		}else if(treestring[i]==')'){
			end = i;
			return (string)"bracket";
		}

	}

	return "error";
}		/* -----  end of function whats_right  ----- */



/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  whats_left
 *  Description:  what is left of the name in the Newick string
 * =====================================================================================
 */
string whats_left (int& start, string tree){

	for(int i=start;i>=0;--i){
		if(tree[i]==','){
			start = i;
			return (string)"comma";
		}else if(tree[i]=='('){
			start = i;
			return (string)"bracket";
		}
	}
	return "error";
}		/* -----  end of function whats_left  ----- */




/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  Put_distances_on_tree
 *  Description:  
 * =====================================================================================
 */


vector<Node *> Put_distances_on_tree(TreeTemplate<Node> *tree, std::auto_ptr<DistanceMatrix> DS ){	

	vector<int> nodes = tree->getNodesId();
	tree->setBranchLengths(0.0);

	//vector<string> nams;
	//	nams = DS->getNames();
	map<string, int> present;
	vector<int> tips;
	vector<Node *> pres;
	vector<int> fathers;
	vector<Node *> fath;

	vector<string> names = tree->getLeavesNames();


	/* make nested nodes using the vector of ints of leaves */
	vector<int> isLeaf;
	isLeaf = tree->getLeavesId();





	/*put the distances on the tips */
	vector<int> temp_distances;
	for(unsigned int j=0;j<isLeaf.size()-2;++j){
		pres.push_back(tree->getNode(isLeaf[j]));
		double dist1, dist2, dist3, dist4;

		dist1 = (*DS)(tree->getNodeName(isLeaf[j]), tree->getNodeName(isLeaf[j+1]));
		dist2 = (*DS)(tree->getNodeName(isLeaf[j]), tree->getNodeName(isLeaf[j+2]));
		dist3 = (*DS)(tree->getNodeName(isLeaf[j+1]), tree->getNodeName(isLeaf[j+2]));
		dist4 = (dist1+dist2-dist3)/2;
		tree->setDistanceToFather(isLeaf[j], dist4);
		temp_distances.push_back(dist4);
		dist4 = (dist1-dist2+dist3)/2;
		tree->setDistanceToFather(isLeaf[j+1], dist4);
		temp_distances.push_back(dist4);
		dist4 = (-dist1+dist2+dist3)/2;
		tree->setDistanceToFather(isLeaf[j+2], dist4);
		temp_distances.push_back(dist4);
	}
	pres.push_back(tree->getNode(isLeaf[isLeaf.size()-2]));
	pres.push_back(tree->getNode(isLeaf[isLeaf.size()-1]));

	vector<int> inner;
	vector<Node *> innod;

	innod = tree->getInnerNodes();
	inner = TreeTemplateTools::getInnerNodesId(*tree->getRootNode());

	map<int, int> idmap;
	map<Node *, int> nodmap;
	int j=0;
	for(unsigned int i=0;i<inner.size();++i){
		if(inner[i]!=tree->getRootId()){

			int da = tree->getFatherId(inner[i]);
			vector<int> sol = tree->getSonsId(da);
			if(sol.size()>1){
				idmap[inner[i]] = j;
				++j;
			}else{
				tree->setDistanceToFather(sol[0], 0.0);
			}
		}
	}

	vector<int> sonny = tree->getSonsId(tree->getRootId());
	vector< vector<double> > pol;
	vector< double> temper(idmap.size(), 0.0);
	vector< double> ans(idmap.size(), 0.0);
	for(unsigned int i=0;i<idmap.size();++i){
		pol.push_back(temper);
	}


	map<int, int> combo;
	int k=0;

	map<int, int>::iterator mi;

	for(mi=idmap.begin();mi!=idmap.end();++mi){

		get_unique_pair(combo, mi->first, tree, idmap, pol, k, ans, DS);
		++k;
	}

	gsl_matrix *A; 
	gsl_vector *x, *b;

	A= gsl_matrix_alloc(pol.size(), pol[0].size());
	b= gsl_vector_alloc(pol[0].size());
	x= gsl_vector_alloc(pol[0].size());


	for(unsigned int i=0;i<pol.size();++i){
		for(unsigned int j=0;j<pol[i].size();++j){
			gsl_matrix_set(A, i, j, pol[i][j]);
		}
		gsl_vector_set(b, i, ans[i]);
	}

	gsl_linalg_HH_solve (A, b, x);

	for(unsigned int i=0;i<pol.size();++i){

		if(gsl_vector_get(x, i)<=0){
			gsl_vector_set(x, i, 0.01);
		}


		if((double)gsl_vector_get(x,i)<=0.0){
			tree->setDistanceToFather(inner[i], gsl_vector_get(x, i));	
		}else{
			tree->setDistanceToFather(inner[i], 0.001);	
		}

	}


	string tempp = TreeTemplateTools::treeToParenthesis(*tree, false);




	for(unsigned int i=0;i<nodes.size();++i){

		if(tree->hasDistanceToFather(nodes[i])==false){
			if(nodes[i]==tree->getRootId()){
				tree->setDistanceToFather(nodes[i], 0.0);
			}else{
				cerr << "Error setting distances on tree!" << endl;
				exit(-1);
			}
		}

	}

	gsl_vector_free(b);
	gsl_vector_free(x);
	gsl_matrix_free(A);


	return pres;
}
/* -----  end of function Put_distances_on_tree  ----- */



/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  get_unique_pair
 *  Description:  
 * =====================================================================================
 */
void get_unique_pair (map<int, int>& combo, int node, TreeTemplate<Node>* tree, map<int, int>& idmap, vector< vector<double> >& pol, int l, vector<double>& ans, std::auto_ptr<DistanceMatrix> DS ){

	Node *subtreeRoot = TreeTemplateTools::cloneSubtree<Node>(*tree->getNode(node));
	TreeTemplate<Node> subtree(subtreeRoot); //Create a new (independent) tree from the subtree
	vector<int> leavesid = subtree.getLeavesId();


	map<int, int>::iterator myit;
	myit = combo.find(leavesid[0]);

	if(myit!=combo.end() && myit->second!=leavesid[leavesid.size()-1]){
		combo[myit->first]=leavesid[leavesid.size()-1];
	}

	int g=0;

	if(tree->hasNodeName(leavesid[0])==false && leavesid.size()==1){
		pol.erase(pol.begin()+l, pol.begin()+l+1);
		return;
	}

	while(tree->hasNodeName(leavesid[g])!=true){
		++g;
	}


	int h = leavesid.size()-1;

	while(tree->hasNodeName(leavesid[h])!=true){
		--h;
		if(h==0){
			return;
		}
	}

	int dar1, dar2;

	dar1=tree->getFatherId(leavesid[g]);
	map<int, int>::iterator mit;
	mit = idmap.find(dar1);
	if(dar1!=node)
		pol[l][mit->second]=1.0;
	while(dar1!=node){
		dar1=tree->getFatherId(dar1);
		if(dar1!=node)
			pol[l][mit->second]=1.0;
	}



	dar2=tree->getFatherId(leavesid[h]);
	mit = idmap.find(dar2);
	if(dar2!=node)
		pol[l][mit->second]=1.0;
	while(dar2!=node){
		dar2=tree->getFatherId(dar2);
		if(dar1!=node)
			pol[l][mit->second]=1.0;
	}


	pol[l][idmap[node]]=1.0;
	ans[l] -= tree->getDistanceToFather(leavesid[g]);
	ans[l] -= tree->getDistanceToFather(leavesid[h]);


	ans[l] += (*DS)(tree->getNodeName(leavesid[g]), tree->getNodeName(leavesid[h]));
}		/* -----  end of function get_unique_pair  ----- */


