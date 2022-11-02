FROM ubuntu:22.04

# Initialise apt-get
RUN apt-get update --fix-missing

# Install programs needed for the server to run
RUN apt-get install -y net-tools
RUN apt-get install -y wget 
RUN apt-get install -y sudo
RUN apt-get install -y vim-tiny
RUN apt-get install -y libgl1-mesa-glx
RUN apt-get install -y cron
RUN apt-get install -y clustalo
RUN apt-get install -y  mafft
RUN apt-get install -y python3-pip
RUN pip install panel bokeh
RUN pip install nglview
RUN pip install biotite
#RUN pip install pybioviz
RUN pip install plotly==5.11.0
RUN pip install dash

# Set the working directory
RUN mkdir -p /opt/MSALigMap
WORKDIR /opt/MSALigMap

# Install XAMPP
#RUN wget https://downloadsapachefriends.global.ssl.fastly.net/8.1.6/xampp-linux-x64-8.1.6-0-installer.run
RUN wget https://sourceforge.net/projects/xampp/files/XAMPP%20Linux/8.1.6/xampp-linux-x64-8.1.6-0-installer.run
RUN chmod +x xampp-linux-x64-8.1.6-0-installer.run
RUN sudo /opt/MSALigMap/xampp-linux-x64-8.1.6-0-installer.run

# Install Anaconda
RUN wget https://repo.anaconda.com/archive/Anaconda3-2022.05-Linux-x86_64.sh
RUN chmod +x Anaconda3-2022.05-Linux-x86_64.sh
RUN ./Anaconda3-2022.05-Linux-x86_64.sh -b -p /usr/local/Anaconda3

# Install python dependencies
RUN /usr/local/Anaconda3/bin/conda install -y -c conda-forge biopython

# Configure a con job to regularly remove SVG files
RUN echo "* 0 * * * root rm /opt/lampp/htdocs/MSALigMap/tmp/*.svg" >> /etc/crontab
RUN echo "* 0 * * * root rm /opt/lampp/htdocs/MSALigMap/tmp/*.zip" >> /etc/crontab

# Install MSALigMap
RUN mkdir /opt/lampp/htdocs/MSALigMap
RUN mkdir /opt/lampp/htdocs/MSALigMap/tmp
#RUN mkdir /opt/lampp/htdocs/MSALigMap/jalview-js
RUN chmod 775 /opt/lampp/htdocs/MSALigMap/tmp
#RUN chmod 775 /opt/lampp/htdocs/MSALigMap/jalview-js
RUN chmod 775 /opt/lampp/htdocs/MSALigMap
RUN chgrp -R daemon /opt/lampp/htdocs/MSALigMap

WORKDIR /opt/lampp/htdocs/MSALigMap
ADD httpd.conf /opt/lampp/etc/
ADD *.py /opt/lampp/htdocs/MSALigMap/
#ADD *.gif /opt/lampp/htdocs/MSALigMap/
#ADD *.png /opt/lampp/htdocs/MSALigMap/
ADD *.jpg /opt/lampp/htdocs/MSALigMap/
ADD *.jpeg /opt/lampp/htdocs/MSALigMap/
#ADD jalview-js /opt/lampp/htdocs/MSALigMap/jalview-js
RUN chmod +x /opt/lampp/htdocs/MSALigMap/*.py
#RUN ln -s /usr/bin/clustalo 
#RUN pip install biotite
#RUN sudo  apt-get install -y  mafft
#RUN sudo apt-get install -y dssp


# TO BE MOVED
#RUN apt-get install -y nano

# Start the server
ENV PATH=/usr/local/Anaconda3/bin:/opt/lampp/htdocs/MSALigMap:$PATH
EXPOSE 80 
CMD /opt/lampp/xampp start; service cron start; bash 
