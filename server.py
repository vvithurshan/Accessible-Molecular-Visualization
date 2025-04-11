from pymol import pymolhttpd

## Pymol httpd server
httpd = pymolhttpd.PymolHttpd(8080, "htdocs")
httpd.start()

## Function to be exposed
httpd.expose("Welcome", Welcome)
httpd.expose("load", load)
httpd.expose("delete", delete)
httpd.expose("zoomin", zoomin)
httpd.expose("zoomout", zoomout)
httpd.expose("rotate_x", rotate_x)
httpd.expose("rotate_y", rotate_y)
httpd.expose("rotate_z", rotate_z)
httpd.expose("reset", reset)
httpd.expose("orient", orient)
httpd.expose("sphere", sphere)