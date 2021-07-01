## Project Information

Team number: xohw21-230  	<br />
Project name: WAPA		<br />
Date: 30/06/2021			<br />
Version of uploaded archive:1	<br />
													<br />
University name: Politecnico di Milano				<br />
Supervisor name: Marco Domenico Santambrogio		<br />
Supervisor e-mail: marco.santambrogio@polimi.it		<br />
Participant: Sofia Breschi							<br />	
Email: sofia.breschi@mail.polimi.it					<br />
Participant: Leonardo De Grandis						<br />
Email: leonardo1.degrandis@mail.polimi.it					<br />
Participant: Alessia Della Valle					<br />
Email: alessia.dellavalle@mail.polimi.it			<br />
													<br />
Board used: Xilinx Alveo U280     <br />
Vitis Version: 2020.2							<br />

## Project Description
<p align="justify">
Genomics is a subdiscipline of biology dedicated to analyzing the functionsof genomes and developing state of art methods to map and sequence it.In the 1970s, recombinant DNA, gene cloning, and DNA sequencing tech-nologies began to flourish. That is why, in the ‘80s, scientists began a dis-cussion about the possibility of sequencing the entire human genome andthis led to the launching of the Human Genome Project (HGP). HGP aimedto sequence the human DNA, mapping human genes, and creating a mapof the human genome. In 1999, the first sequenced human chromosome waspublished; in 2004 an almost complete sequence of the human genome wasreleased.  Since then, rapid progress was made in this emerging field,making the use of genomic data more and more common. In the last fewyears, the DNA-sequencing procedure cost is characterized by a decreasingtrend. This tendency can lead to increased utilization in medicine of NextGeneration Sequencing (NGS), in particular for the diagnosis of diseases andthe production of personalized drugs.
For a vast field of applications, looking for a pattern in a reference sequencelike the entire human genome, may result in a great achievement for devel-oping further studies.Pattern matching problems aim to find the occurrences of a given se-quence into a reference one.  This can be accomplished with various ap-proaches, from the exact matching methods, to gap-affine methods whichevaluate the results based on a scoring system.  Exact matching providesa better result but at the cost of being much more computationally inten-sive and time-consuming. Furthermore, these methods are not optimized forgeneral purpose architectures because of the huge amount of data requiredwhen processing.Our work consist of an implementation on Field Programmable GateArray of the BWA-MEM algorithm, a widely-used tool for genomic shortread mapping.  The aim of the implementation is to improve the overallperformance of one of the most time-consuming stages of the alignment,which is based on a banded gap affine Smith–Waterman algorithm.The application will run on a Xilinx Alveo U280 board.
</p>

## Usage

The repo already includes a host and a bitstream for the Alveo U280.
Before executing WAPA be sure so source XRT.

### Compilation
Before executing WAPA be sure so source XRT and Vitis.
WAPA requires Vitis 2020.2. To build WAPA simply type:
```
make all TARGET=hw
```
WAPA has been written to run on the Xilinx Alveo U280.
WAPA will an executable called **host** and a bitstream file for the Alveo U280 called **SW.xclbin**.

### Demo

To check everything works properly type:
```
./host build_dir.hw.xilinx_u280_xdma_201920_3/SW.xclbin 
```
This command executes WAPA on a randomly generated datasets with sequences of 1024 read pairs.
