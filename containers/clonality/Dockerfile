FROM python:3.11-slim
LABEL org.opencontainers.image.source="https://github.com/trethewey/lymphix"
LABEL org.opencontainers.image.description="Clonality metrics + HTML report generator"

RUN pip install --no-cache-dir \
        pandas==2.2.2 \
        numpy==1.26.4 \
        scipy==1.13.1 \
        jinja2==3.1.4 \
        plotly==5.22.0

COPY bin/clonality_metrics.py /opt/scripts/clonality_metrics.py
COPY bin/generate_report.py   /opt/scripts/generate_report.py
COPY bin/report_template.html /opt/scripts/report_template.html
RUN chmod +x /opt/scripts/*.py
ENV PATH=/opt/scripts:$PATH

WORKDIR /work
