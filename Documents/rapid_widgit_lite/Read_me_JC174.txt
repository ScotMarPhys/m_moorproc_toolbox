Note that mooring positions are in moor_pos.dat
This is a symbolic link normally link to moor_pos_target.dat,
but some times change. E.g. used when tracking a float during JC174

There is a version rapid_usbl_lite which also takes data in from teh USBL data stream
for use when draggin or deploying MYRTLE etc

Can now set a fallback distance (m) as an argument after the mooring name.
For this to work correctly assumes that the ship is heading towards the mooring.

Now runnning on workstations Koaekea CnetOS 7
To start and check webserver use
- systemctl start httpd
- systemctl status httpd
May need ot open firewall
 - firewall-cmd --zone=public --add-service=http
However, shoudl autommatically start as have used
 - systemctl enable httpd
 - firewall-cmd  --permanent  --zone=public --add-service=http

See previous notes.