import socket
from grex_t2 import socket_grex
import logging as logger

logger.basicConfig(filename="output.log", encoding="utf-8", level=logger.DEBUG)

HOST = "127.0.0.1"
PORT = 12345

def main(trigger=True):
    # Use roughly 8 seconds as a gulp size
    gulpsize = 16384 * 8

    s = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
    # Create a UDP socket
    s = socket.socket(socket.AF_INET, socket.SOCK_DGRAM)
    # Ensure that you can reconnect
    s.setsockopt(socket.SOL_SOCKET, socket.SO_REUSEADDR, 1)
    # Bind the socket to the port
    server_address = (HOST, PORT)
    s.bind(server_address)

    print("Connected to socket %s:%d. Triggering set to %s" % (HOST, PORT, trigger))
    logger.info("Connected to socket %s:%d. Triggering set to %s" % (HOST, PORT, trigger))

    candsfile = ["", "", "", "", ""]

    # Outer loop that runs as long as T2 is running
    while True:
        candstr_list = ''
        cand_count = 0
        # Inner loop for chunks of Heimdall output data
        while True:
            # Recieve 512 bytes
            data, address = s.recvfrom(512)

            # Waiting for end of text. When chunk is done, break inner loop
            if len(data)==1 and data==b'\x03':
                break
            
            candstr = data.decode("utf-8")
            
            # Removing the \n from the end of line
            candstr_list += candstr
            cand_count += 1

        print("Number of candidates %d" % cand_count)
        print(candstr_list)
        if cand_count > 0:
            print("Filtering")
            socket_grex.filter_candidates(candstr_list, trigger=trigger)
            print("Finished filtering")
            
        continue
        
        try:
            itime = int(candstr.split("\t")[2])
        except IndexError:
            continue
            
        # Read time sample to keep track of gulp number
        itime = int(candstr.split("\t")[2])
        gulp_ii = itime // gulpsize

        if candsfile == ["", "", "", "", ""]:
            gulp = gulp_ii
            print("Starting gulp is %d" % gulp)
            logger.info("Starting gulp is %d" % gulp)
        if gulp_ii - gulp < 0:
            print("Receiving candidates gulps from before current gulp")
            logger.info("Receiving candidates gulps from before current gulp")            
            print(gulp_ii, gulp)
            continue
        if gulp_ii - gulp >= len(candsfile):
            print("Receiving candidates too far ahead of current gulp")
            logger.info("Receiving candidates too far ahead of current gulp")
            print(gulp_ii, gulp)
            continue
        candsfile[gulp_ii - gulp] += candstr

        # If cand is received with gulp 2 or more than
        # current gulp, process current gulp
        if gulp_ii >= gulp + 3:
            if candsfile[0] == "":
                candsfile.pop(0)
                candsfile.append("")
                continue
            gulp += 1
            socket_grex.filter_candidates(candsfile[0], trigger=trigger)
            candsfile.pop(0)
            candsfile.append("")
            continue

    exit()
