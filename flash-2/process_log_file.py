#!/usr/bin/env python

def process_log_file(logfile):
    import datetime as dt

    logstamps = []
    starttimes = []
    stoptimes = []
    ncpus = []
    nsteps = 0

    # Open the log file
    f = open(logfile,'r')

    # Read the first line
    line = f.readline()

    # Process the opening information
    runinit = line.strip().split()           # declaration of new sim start
    sd = runinit[3].split('-')               # split the date 'dd-mm-yyyy'
    sd = [int(x) for x in sd]                # convert to integers
    st = runinit[4].split(':')               # split the time 'hh:mm.mm'
    st.append('')                            # extend the start time array
    st[1],st[2] = divmod(float(st[1])*60,60) # convert the fractional minute to minutes and seconds
    st = [int(x) for x in st]                # convert to integers
    startdt = dt.datetime(sd[2],sd[0],sd[1],st[0],st[1],st[2]) # create datetime object
    logstamps.append(startdt)               # add to log stamps array

    line = f.readline()
    while line != '':
        # Capture the current position and advance the file one line 
        last_pos = f.tell()
        line = f.readline()
        if line == '':
            # The party is over, we've reached the end of the file
            evostop = lastline.split()
            sd = evostop[1]
            st = evostop[2]
            sd = sd.split('-')                       # split the date 'dd-mm-yyyy'
            sd = [int(x) for x in sd]                # convert to integers
            st = st.split(':')                       # split the time 'hh:mm.mm'
            st.append('')                            # extend the start time array
            st[1],st[2] = divmod(float(st[1])*60,60) # convert the fractional minute to minutes and seconds
            st = [int(x) for x in st]                # convert to integers
            stopdt = dt.datetime(sd[2],sd[0],sd[1],st[0],st[1],st[2]) # create datetime object
            stoptimes.append(stopdt)
        if 'FLASH log file' in line:
            # Since this shouldn't be the first initialization of a FLASH run,
            # rewind the log file a little and capture the datetime stamp to
            # figure out how long it ran
            f.seek(last_pos)
            evostop = lastline.split()
            sd = evostop[1]
            st = evostop[2]
            sd = sd.split('-')                       # split the date 'dd-mm-yyyy'
            sd = [int(x) for x in sd]                # convert to integers
            st = st.split(':')                       # split the time 'hh:mm.mm'
            st.append('')                            # extend the start time array
            st[1],st[2] = divmod(float(st[1])*60,60) # convert the fractional minute to minutes and seconds
            st = [int(x) for x in st]                # convert to integers
            stopdt = dt.datetime(sd[2],sd[0],sd[1],st[0],st[1],st[2]) # create datetime object
            stoptimes.append(stopdt)
            f.readline()
            
            # Now handle the new initialization
            runinit = line.strip().split()           # declaration of new sim start
            sd = runinit[3].split('-')               # split the date 'dd-mm-yyyy'
            sd = [int(x) for x in sd]                # convert to integers
            st = runinit[4].split(':')               # split the time 'hh:mm.mm'
            st.append('')                            # extend the start time array
            st[1],st[2] = divmod(float(st[1])*60,60) # convert the fractional minute to minutes and seconds
            st = [int(x) for x in st]                # convert to integers
            startdt = dt.datetime(sd[2],sd[0],sd[1],st[0],st[1],st[2]) # create datetime object
            logstamps.append(startdt)               # add to start times array
        if 'Number of processors' in line:
            ncpus.append(int(line.split(':')[1].strip()))
        if 'Build stamp' in line:
            buildtime = line.split(':')[1:]
            buildtime = ''.join(buildtime).strip()
        if 'System info' in line:
            sysinfo = line.split(':')[1].strip()
        if 'Enter evolution loop' in line:
            evostart = line.split()
            sd = evostart[1]
            st = evostart[2]
            sd = sd.split('-')                       # split the date 'dd-mm-yyyy'
            sd = [int(x) for x in sd]                # convert to integers
            st = st.split(':')                       # split the time 'hh:mm.mm'
            st.append('')                            # extend the start time array
            st[1],st[2] = divmod(float(st[1])*60,60) # convert the fractional minute to minutes and seconds
            st = [int(x) for x in st]                # convert to integers
            startdt = dt.datetime(sd[2],sd[0],sd[1],st[0],st[1],st[2]) # create datetime object
            starttimes.append(startdt)
        if 'step:' in line:
            stepline = line.split()
            nsteps = int(stepline[5][2:])
        lastline = line

    # Close the file
    f.close()

    # Calculate runtimes and total cpuhours
    runtimes = []
    cpuhours = 0.0
    wallclock = 0.0
    for i in range(len(stoptimes)):
        a = stoptimes[i]-starttimes[i]
        hours = a.days*24.0+a.seconds/60.0/60.0
        runtimes.append(hours)
        wallclock += runtimes[i]
        cpuhours += runtimes[i]*ncpus[i]

    loginfo = {}
    loginfo['buildtime'] = buildtime
    loginfo['sysinfo'] = sysinfo
    loginfo['runtimes'] = runtimes
    loginfo['wallclock'] = wallclock
    loginfo['cpuhours'] = cpuhours
    loginfo['nsteps'] = nsteps

    return loginfo
