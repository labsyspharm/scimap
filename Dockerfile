FROM python:3.8

RUN pip install scimap

COPY . /app/