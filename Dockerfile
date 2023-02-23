FROM python:3.8

RUN pip install --no-cache-dir scimap --upgrade

COPY . /app/