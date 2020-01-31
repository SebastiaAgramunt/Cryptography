# Tutorial on Cryptogarphy and Multiparty Computation

This is a quick introduction to cryptography based on the books [An introduction to mathematical cryptograpy](https://www.springer.com/gp/book/9781441926746) by Hoffstein, Pipher and Silverman and [Introduction to modern cryptography](https://www.crcpress.com/Introduction-to-Modern-Cryptography/Katz-Lindell/p/book/9781466570269) by Katz and [Lindell](https://u.cs.biu.ac.il/~lindell/). 

The main objective of this repository is to show from first principles some of the mathematical foundations underlying modern cryptography. At the same time implementing some algorithms to achieve a certain degree of security. 

Once some concepts on cryptography are clear we show more complex concepts like mulitparty computation.

## Getting Started

These instructions will get you a copy of the project up and running on your local machine. 

### Prerequisites

Make sure you have installed [Docker](https://www.docker.com/get-started) in your computer. Try to get your docker version on the command line

```sh 
docker --version
```
It's been tested on  ```Docker version 19.03.5``` but upper versions may work as well.


### Running

In the main folder run

```sh
make build-run
```

This will create the image and run the container. By default image name is *cryptography* and container name *crypt* (you can change this in the Makefile).

## Accessing the notebooks

After making the *build-run* you should have a docker container named *crypt* running. To access the notebooks just open in your browser

```
http://localhost:8888/
```

you can stop, run and delete the image by typing

```
make stop
make run
make delete
```

delete command will erase both, the docker image and the container.