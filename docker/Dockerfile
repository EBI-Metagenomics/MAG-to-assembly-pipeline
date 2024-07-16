FROM python:3.9-slim

RUN pip install --no-cache-dir argparse requests xmltodict tqdm retry biopython

COPY bin/ /usr/local/bin/

WORKDIR /usr/local/bin

CMD ["python3"]