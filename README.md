AutoCoEV

# AutoCoEv – a high-throughput in silico pipeline for predicting inter-protein co-evolution

_Petar B. Petrov<sup>1,2*</sup>, Luqman O. Awoniyi<sup>1,2</sup>, Vid Šuštar1<sup>1</sup>, M. Özge Balcı<sup>1,2</sup> and Pieta K. Mattila<sup>1,2</sup>_

1. Institute of Biomedicine and MediCity Research Laboratories, University of Turku, Turku, Finland
2. Turku Bioscience, University of Turku and Åbo Akademi University, Turku, Finland

**Abstract**

Protein-protein communications govern cellular processes via complex regulatory networks but our understanding of them represents merely a tip of an iceberg. Thus, identifying novel interactions between proteins can significantly facilitate our comprehension of the mechanistic principles of protein functions. Co-evolution between proteins is a sign of functional communication and, as such, provides a powerful approach to search for novel direct or indirect molecular partners. However, evolutionary analysis of large arrays of proteins, in silico, is a highly time-consuming effort, which has limited the usage of this method to protein pairs or small protein groups. Here, we developed AutoCoEv, a user-friendly computational pipeline for the search of co-evolution between a large number of proteins. By driving 15 individual programs, culminating in CAPS2 as the software for detecting co-evolution, AutoCoEv achieves seamless automation and parallelization of the workflow. Importantly, we provide a patch to CAPS2 source code to strengthen its statistical output, allowing for multiple comparisons correction and enhanced analysis of the results. We apply the pipeline to inspect co-evolution among 324 individual proteins identified to locate at the lipid rafts of B lymphocytes. We successfully detected multiple strong coevolutionary connection between the proteins, predicting many novel partners and previously unidentified clusters of functionally related molecules. 

--
* Requests for the scripts should be addressed to Petar Petrov (pebope ат utu.fi)


**The 'main' branch is currently under development and things may break.**
**Get the latest release if you want to give AutoCoEv an try!**
