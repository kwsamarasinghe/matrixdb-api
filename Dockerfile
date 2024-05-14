FROM python:3.11
WORKDIR /usr/src/app

COPY requirements.txt ./
RUN pip3.11 install --no-cache-dir -r requirements.txt

COPY . .

ENV PYTHONPATH="${PYTHONPATH}:./src"
CMD ["uwsgi", "--ini", "uwsgi.ini", "--http", ":8000", "--wsgi-file", "src/matrixdb/api/main.py", "--callable", "app", "--processes", "2", "--threads", "2"]