FROM python:3.11-slim
LABEL org.opencontainers.image.source="https://github.com/trethewey/lymphix"
LABEL org.opencontainers.image.description="Clonality metrics + HTML report generator"

RUN pip install --no-cache-dir \
        pandas==2.2.2 \
        numpy==1.26.4 \
        scipy==1.13.1 \
        jinja2==3.1.4 \
        plotly==5.22.0

# The whole of bin/, not an enumerated list. The scripts import each other by
# bare module name (Nextflow puts <projectDir>/bin on the task PATH, and the
# shared definitions live in bin/lymphix_common.py), so naming files
# individually means a new module is simply absent inside the image and the
# ImportError surfaces only when a real sample reaches the task. The list here
# was already wrong: it copied bin/report_template.html, which does not exist,
# so this image could not be built at all.
COPY bin/ /opt/scripts/

# The reports inline the logo as SVG. Without the assets the loader falls back
# to a text wordmark, so container-generated reports quietly lost their
# branding while locally generated ones kept it. Only the SVGs are needed;
# lymphix_common.load_logo_svg() looks in <script dir>/../assets.
COPY assets/*.svg /opt/assets/

RUN chmod +x /opt/scripts/*.py
ENV PATH=/opt/scripts:$PATH

WORKDIR /work
