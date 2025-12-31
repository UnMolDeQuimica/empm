FROM python:3.12-slim

# Port used by this container to serve HTTP.
EXPOSE 8005

# Set environment variables
ENV PYTHONUNBUFFERED=1 \
    PORT=8006

# Install system dependencies for RDKit rendering
RUN apt-get update && apt-get install -y \
    libxrender1 \
    libxext6 \
    libexpat1 \
    libsm6 \
    libglib2.0-0 \
    libgomp1

# Set work directory
WORKDIR /code

# Install dependencies
COPY requirements.txt /code/
RUN pip install --upgrade pip
RUN pip install -r requirements.txt

# Copy project
COPY . /code/

# Run migrations and start server
CMD python manage.py migrate && python manage.py runserver 0.0.0.0:8006
