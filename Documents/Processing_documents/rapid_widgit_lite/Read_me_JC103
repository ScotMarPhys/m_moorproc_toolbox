How to display navigation data (and other data) produced by Matlab on a webserver while at sea

DAS 16 October 2012 on D382
updated
DAS April on JC10

0)  information about mooring locations is in the file moor_pos.dat
I creaated three versions. Control which one is used by symbolic link
 - moor_pos_target.dat 	- Target positions for mooring deployments.  This is one used during a deployment.
 - moor_pos_jc103_dep.dat  - actual positions of deployed moorigns after triangulation
 - moor_pos_jc103_rec.dat  - actual posotions of moorings to be recovered
 
1)  To create the plots
The script now renamed rapid_widgit_lite.m requires one arguement which is the name of the mooring. 
If use 'ship' instead of a mooring a plot based on current position will be made.
There are also optional arguements that are nto normally needed. These are for size of plot 
and contour interval- see the help lines in the file.
Progam calls one other routine - rapid_web_map.m.
There is a loop in the program and a plot is created every time interval and saved in this
directory to: latest_plot.png.

Check that the script has the correct positions!!

Should only be run in one matlab session at any time, otherwise plot files will be overwritten.

2) Start matlab   
The command figure('visible', off') is used so script will run without plotting to the screen.
It is best to run in a terminal with -nodesktop option for matlab.   

3) Setup a webserver
On this  cruise changed so that webserver is now on workstation Banba.
initially PHP did not work and found that there was a missing Apache module. After finding version of PHP with
> rpm -q php53
Installed module with
>rpm -ivh apache2-mod_php53-5.3.17-0.13.7.x86_64.rpm

Used YaST, System -> System services  and enable run level 3, 5 so shoudl start up on boot.  
Manually turn on off with 
/usr/sbin/rcapache2 start (stop)

Created a virtual server by adding a .conf file in the directory:
/etc/apache2/vhosts.d
Based on the template file and made changes so that
   ServerAdmin das@noc.ac.uk
   ServerName banba.local
   DocumentRoot /srv/www/vhosts/jc103
Latter is where must place files

) Create a webpage
Place files in: /srv/www/vhosts/jc103

Create a file index.php. The new file shoudl include lines

<?php header('Refresh: 10'); ?>
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

The first line will cause the  page to reload every 10 seconds in a browser.
So last step is to get a copy of the matlab plot in this directory.

7)  View the webpage
On JC103 the URL was:
http://banba.local
