# NGSadmix
Inferring admixture proportions from NGS data

see website
http://www.popgen.dk/software/index.php/NgsAdmix

# Tutorial

For a repository-local tutorial using the bundled demo input files, see [TUTORIAL.md](TUTORIAL.md).

# INSTALL

To install clone and compile. Only tested using GNU compiler

## clone

```
git clone https://github.com/aalbrechtsen/NGSadmix.git
```

## compile
```
cd NGSadmix
g++ NGSadmix.cpp -O3 -lpthread -lz -o NGSadmix
```

