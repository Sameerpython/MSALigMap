FROM ubuntu:20.04

# Initialise apt-get
RUN apt-get update

# Install programs needed for the server to run
RUN apt-get install -y net-tools
RUN apt-get install -y wget 
RUN apt-get install -y sudo
RUN apt-get install -y vim-tiny
RUN apt-get install -y libgl1-mesa-glx
RUN apt-get install -y cron
RUN apt-get install -y clustalo
RUN apt-get install -y python3.9


RUN mkdir -p /opt/MSALigMap
WORKDIR /opt/MSALigMap

# Install XAMPP
RUN wget https://www.apachefriends.org/xampp-files/7.3.9/xampp-linux-x64-7.3.9-0-installer.run
RUN chmod +x xampp-linux-x64-7.3.9-0-installer.run
RUN sudo /opt/bindingdata/xampp-linux-x64-7.3.9-0-installer.run

# Install Anaconda
#RUN wget https://repo.anaconda.com/archive/Anaconda2-2019.10-Linux-x86_64.sh
#RUN chmod +x Anaconda2-2019.10-Linux-x86_64.sh
#RUN sudo ./Anaconda2-2019.10-Linux-x86_64.sh -b -p /usr/local/Anaconda2.7
#RUN echo "PATH=/usr/local/Anaconda2.7/bin:$PATH" >> /etc/profile

# Install python dependencies
#RUN /usr/local/Anaconda2.7/bin/conda install -y -c conda-forge biopython
#RUN export PATH=/usr/local/Anaconda2.7/bin:$PATH

# Configure a con job to regularly remove SVG files
RUN echo "* 0 * * * root rm /opt/lampp/htdocs/MSALigMap/tmp/*.svg" >> /etc/crontab
RUN echo "* 0 * * * root rm /opt/lampp/htdocs/MSALigMap/tmp/*.zip" >> /etc/crontab

# Install MSALigMap
RUN mkdir /opt/lampp/htdocs/MSALigMap
RUN mkdir /opt/lampp/htdocs/MSALigMap/tmp
RUN chmod 775 /opt/lampp/htdocs/MSALigMap/tmp
RUN chmod 775 /opt/lampp/htdocs/MSALigMap
RUN chgrp -R daemon /opt/lampp/htdocs/MSALigMap

WORKDIR /opt/lampp/htdocs/MSALigMap
ADD httpd.conf /opt/lampp/etc/
ADD *.py /opt/lampp/htdocs/MSALigMap/
#ADD *.gif /opt/lampp/htdocs/MSALigMap/
#ADD *.png /opt/lampp/htdocs/MSALigMap/
#ADD *.jpg /opt/lampp/htdocs/MSALigMap/
#ADD *.jpeg /opt/lampp/htdocs/MSALigMap/
RUN chmod +x /opt/lampp/htdocs/MSALigMap/*.py

# Start the server
EXPOSE 80 
CMD /opt/lampp/xampp start; service cron start; bash 
