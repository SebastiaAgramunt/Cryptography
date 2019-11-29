FROM python:3.7.5
#FROM python:3.7.2-stretch

# Updating repository sources
RUN apt-get update

COPY requirements.txt /tmp/
RUN pip install --upgrade pip
RUN pip install --requirement /tmp/requirements.txt

# Installing requirements
RUN apt-get install cron -yqq \
    curl

# Setting Working Directory
WORKDIR /home

# Launching Jupyter Notebook
# change token for more security
EXPOSE 8888
CMD jupyter notebook --no-browser --ip=0.0.0.0 --allow-root --NotebookApp.token='' --NotebookApp.password=''
