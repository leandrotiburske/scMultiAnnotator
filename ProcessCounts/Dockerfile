FROM python:3.9.20-slim-bullseye

RUN pip install --upgrade pip

RUN pip install --upgrade setuptools

RUN pip install pandas==2.0.0 \
                numpy==1.26.4 \
                'scanpy[leiden]' \
                louvain \
                requests \
                igraph \
                scikit-learn \
                leidenalg \
                argparse

RUN pip install --user scikit-misc
 
CMD [ "/bin/bash" ]