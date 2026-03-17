How to display navigation data (and other data) produced by Matlab on a webserver while at sea

DAS 16 October 2012 on D382

1)  To create the plots
The script rapid_disp_web.m has one arguement which is the name of the mooring. 
If use 'ship' instead of a mooring a plot based on current positino will be made.
This calls one other program das_rapid_map2.m.
There is a loop in the program and a plot is created every one minute and saved in this
directory to: latest_plot.png.

Check that the script has the correct positions!!

Should only be run in one matlab session at any time, otherwise plot files will be overwritten.

2) Start matlab   
The command figure('visible', off') is used so script will run without plotting to the screen.
It is best to run in a terminal with -nodesktop option for matlab.   

3) Setup a webserver
It is easy to set up webserving on a Mac (probably not to hard to do on Linus either).  
In System Preferences see Sharing to turn on Web sharingg.
You will also need to enable PHP.  To do this open a terminal on the Mac and go to
/etc/apache2
and edit the file
httpd.conf
You need to un-comment out the lines:
LoadModule php5_module        libexec/apache2/libphp5.so
LoadModule fastcgi_module     libexec/apache2/mod_fastcgi.so

You will need to restart the Mac to get webserver running with php

4) Create a webpage
Web pages  on a Mac are in 
/Library/WebServer/Documents 

Create a file index.php (remove index.html). The new file shoudl include lines

<?php header('Refresh: 30'); ?>
<!DOCTYPE html PUBLIC "-//W3C//DTD XHTML 1.0 Transitional//EN" "http://www.w3.org/TR/xhtml1/DTD/xhtml1-transitional.dtd">

<html xmlns="http://www.w3.org/1999/xhtml">

        <head>
                <meta http-equiv="content-type" content="text/html;charset=utf-8" />
                <meta name="generator" content="Adobe GoLive" />
                <title>D382 navigation data</title>
                <style type="text/css" media="screen"><!--
                --></style>
        </head>

        <body>
                        <img src="latest_plot.png"  border="0" />
        </body>
</html>

The first line will cause the  page to reload every 30 seconds in a browser.
So last step is get a copy of the matlab plot in this directory.

5) Now you need to get the plotfile onto your webserver
Get Samba running on the Linux workstation and mount on Mac
Use Finder -> Connect to server on the Mac

6) Copy the file to the Webserver  directory.  
I used 'crontab -e ' to create a cron job that runs (on the Mac) every minute and copies the
png file each time.
Don't forget to remove this at the end of the cruise!

7)  View teh webpage
On D382 the URL was:
http://default-dhcp-40.discovery.local
or 
http://default-dhcp-65.discovery.local
