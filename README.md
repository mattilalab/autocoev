# AutoCoEv -- A High-Throughput In Silico Pipeline for Predicting Inter-Protein Coevolution

**NOTE (7 Mar 2022): An updated manual and virtual machines images will be available very soon! Send any requests to Petar Petrov (pebope ат utu.fi)**

## Virtual machine images with everything preinstalled are available here:

* CRUX 3.6.1: [autocoev-crux-3.6.1-20220307.ova](https://seafile.utu.fi/d/a8de85062abf4ab68de9/)
* An OVA VM image with Ubuntu is in preparation

## An updated manual is in preparation (first draft ready!):

* [Manual_AutoCoEv.pdf](doc/Manual_AutoCoEv.pdf)

## Paper

* Published in [Int. J. Mol. Sci. 2022, 23(6), 3351;](https://www.mdpi.com/1422-0067/23/6/3351)

## Abstract
_Petar B. Petrov<sup>1,2*</sup>, Luqman O. Awoniyi<sup>1,2</sup>, Vid Šuštar1<sup>1</sup>, M. Özge Balci<sup>1,2</sup> and Pieta K. Mattila<sup>1,2</sup>_

1. MediCity Research Laboratories, Institute of Biomedicine, University of Turku, 20014 Turku, Finland
2. Turku Bioscience Centre, University of Turku and Åbo Akademi University, 20520 Turku, Finland

Protein–protein interactions govern cellular processes via complex regulatory networks, which are still far from being understood. Thus, identifying and understanding connections between proteins can significantly facilitate our comprehension of the mechanistic principles of protein functions. Coevolution between proteins is a sign of functional communication and, as such, provides a powerful approach to search for novel direct or indirect molecular partners. However, an evolutionary analysis of large arrays of proteins in silico is a highly time-consuming effort that has limited the usage of this method for protein pairs or small protein groups. Here, we developed AutoCoEv, a user-friendly, open source, computational pipeline for the search of coevolution between a large number of proteins. By driving 15 individual programs, culminating in CAPS2 as the software for detecting coevolution, AutoCoEv achieves a seamless automation and parallelization of the workflow. Importantly, we provide a patch to the CAPS2 source code to strengthen its statistical output, allowing for multiple comparison corrections and an enhanced analysis of the results. We apply the pipeline to inspect coevolution among 324 proteins identified to be located at the vicinity of the lipid rafts of B lymphocytes. We successfully detected multiple coevolutionary relations between the proteins, predicting many novel partners and previously unidentified clusters of functionally related molecules. We conclude that AutoCoEv, can be used to predict functional interactions from large datasets in a time- and cost-efficient manner.

