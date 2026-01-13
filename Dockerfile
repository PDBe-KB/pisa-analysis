FROM pdbegroup/pisa AS base

# Copy pisa-analysis code
FROM python:3.9-slim

# Create non-privileged user
RUN useradd -m -u 10001 app \
 && mkdir -p /app && chown -R app:app /app

# Copy pisa-analysis libraray
COPY --chown=app:app pisa_utils /app/pisa_utils
COPY --chown=app:app pisa_utils/run.py /app/run.py
COPY --chown=app:app requirements.txt /app/requirements.txt
COPY --chown=app:app __version__.py /app/__version__.py
COPY --chown=app:app setup.py /app/setup.py
COPY --chown=app:app setup.cfg /app/setup.cfg

WORKDIR /app

# Dependencies
RUN pip install --upgrade pip
RUN pip install -r requirements.txt
RUN pip install .

# Add PISA binary
COPY --from=base /usr/bin/pisa /usr/bin/pisa
ENV PATH=$PATH:/usr/bin/pisa

# Add PISA setup configs
COPY --from=base /usr/share/pisa/setup /usr/share/pisa/setup
ENV PISA_SETUP_DIR=/usr/share/pisa/setup

# Switch to non-privileged user
USER app
