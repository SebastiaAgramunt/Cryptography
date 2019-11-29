.SILENT:
.DEFAULT_GOAL := help

PROJECT := Cryptography (Author Sebastia Agramunt)

IMAGE_NAME=cryptography
CONTAINER_NAME=crypt
PORT=8888

COLOR_RESET = \033[0m
COLOR_COMMAND = \033[36m
COLOR_YELLOW = \033[33m
COLOR_GREEN = \033[32m
COLOR_RED = \033[31m

## Build docker image
build:
	@echo "${COLOR_GREEN}----\nBuilding docker image ${IMAGE_NAME}...\n----\n${COLOR_RESET}"
	docker build -t $(IMAGE_NAME) .

## Run docker image on a container
run:
	@echo "${COLOR_GREEN}----\nBuilding container ${CONTAINER_NAME} and running...\n----\n${COLOR_RESET}"
	docker run -d -v $(shell pwd):/home/ -p 8888:8888 --name $(CONTAINER_NAME) -i $(IMAGE_NAME)

## Build docker image and run it in container
build-run: build run

## Stop docker container
stop:
	@echo "${COLOR_GREEN}----\nStopping container ${CONTAINER_NAME}...\n----\n${COLOR_RESET}"
	docker stop $(CONTAINER_NAME)

## Start the docker container
start:
	@echo "${COLOR_GREEN}----\nStarting container ${CONTAINER_NAME}...\n----\n${COLOR_RESET}"
	docker start $(CONTAINER_NAME)

## Stop docker container and remove containers and images
remove: stop
	@echo "${COLOR_GREEN}----\nStopping container ${CONTAINER_NAME}...\n----\n${COLOR_RESET}"
	docker rm $(CONTAINER_NAME)
	@echo "${COLOR_GREEN}----\nRemoving Image ${IMAGE_NAME}...\n----\n${COLOR_RESET}"
	docker rmi $(IMAGE_NAME)

## Prints help message
help:
	printf "\n${COLOR_YELLOW}${PROJECT}\n------\n${COLOR_RESET}"
	awk '/^[a-zA-Z\-\_0-9\.%]+:/ { \
		helpMessage = match(lastLine, /^## (.*)/); \
		if (helpMessage) { \
			helpCommand = substr($$1, 0, index($$1, ":")); \
			helpMessage = substr(lastLine, RSTART + 3, RLENGTH); \
			printf "${COLOR_COMMAND}$$ make %s${COLOR_RESET} %s\n", helpCommand, helpMessage; \
		} \
	} \
	{ lastLine = $$0 }' $(MAKEFILE_LIST) | sort
	printf "\n"
