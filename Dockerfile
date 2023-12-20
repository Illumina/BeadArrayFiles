FROM python:3.9-slim

WORKDIR /opt/BeadArrayFiles/

# Installing ps to support usage in Nextflow ICA pipelines
RUN apt-get update && apt-get install -y --no-install-recommends procps

COPY . /opt/BeadArrayFiles/

RUN python setup.py install
RUN pip install --no-cache-dir numpy pandas scipy
