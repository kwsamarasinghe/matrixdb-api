FROM python:3.8
WORKDIR /usr/src/app

COPY requirements.txt ./
RUN pip3.8 install --no-cache-dir -r requirements.txt

COPY . .

ENV PYTHONPATH="${PYTHONPATH}:./src"
CMD ["uwsgi", "--http", ":8000", "--wsgi-file", "src/matrixdb/api/main.py", "--callable", "app", "--processes", "20"]