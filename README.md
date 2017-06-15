# clusperin
Tool to cluster proteins (Fasta format).

## Getting Started

**clusperin** is a crucial part of a major set of packages: 

* keggreader
* keggimporter
* anendb


But to run **clusperin** you only need a set of Fasta files with protein sequences and a proper configuration file.

Those Fasta files is expected to be generated by **anendb** web interface.

But you can generate those files manualy. They have to follow the rules below:

* Each file is a set of proteins from a single EC number.
* Example of a file name:

```
EC_2.3.4.5.fasta
```

**Notice** the file format: 'EC\_' **plus** your EC number **plus** '.fasta'

* All the files have to be in a specific directory set in the .clusperin configuration file.
* **ALL** the parameters in the .clusperin.conf file are **mandatory**.
* Example of .clusperin.conf:

```
[clustering]
# Where clusperin expects the Fasta files.
ec_files = /var/kegg/clustering/

# Where clusperin will put result files.
cluster_files = /var/kegg/clustering/clusters/

# The score/cutoff to group proteins.
cutoff = 120 

[log]
# Where to log all the process.
log_file = /var/kegg/clustering/clustering.log
```

* The **.clusprein.conf** file have to be in your $HOME directory.


## Requirements

* BLAST software: Download: ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/. That's the software that actual make the similarity analysis.
* A proper configuration file **.clusperin.conf** in your $HOME directory.


## Installing

* Download (and extract) the zip file from this repository or clone it using **git** command.
* Go to the opened directory.
* Run:

```
python setup.py install
```

* Installation is done.


## Configuring

* Copy the .clusterin.conf example below to your $HOME directory and make the changes that reflects your system.

```
[clustering]
# Where clusperin expects the Fasta files.
ec_files = /var/kegg/clustering/

# Where clusperin will put result files.
cluster_files = /var/kegg/clustering/clusters/

# The score/cutoff to group proteins.
cutoff = 120 

[log]
# Where to log all the process.
log_file = /var/kegg/clustering/clustering.log

```

## Run clusperin

* Simply type **clusperin** in your console.

If your files are correct and all is set properly, **clusperin** will execute all the process without any needs for user's interaction.






