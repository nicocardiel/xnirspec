/*
 * freadline.c
 *
 * Line input through GNU readline (file name completion with TAB and
 * command history) for the Fortran function READC_FILE (readc_file.f90).
 *
 * This routine is free software; you can redistribute it and/or modify it
 * under the terms of the GNU General Public License as published by the Free
 * Software Foundation; either version 3 of the License, or (at your option)
 * any later version. See the file COPYING for details.
 */

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>

#ifdef HAVE_READLINE
#include <readline/readline.h>
#ifdef HAVE_READLINE_HISTORY_H
#include <readline/history.h>
#endif
#endif

/*
 * Read one line from the keyboard using readline.
 *
 * Parameters
 * ----------
 * prompt : null-terminated prompt string, displayed by readline (readline
 *          must know the prompt in order to redraw the line after listing
 *          the possible completions).
 * buffer : Fortran character variable (not null-terminated) that receives
 *          the line, padded with blanks. Lines longer than buflen are
 *          silently truncated, as a Fortran READ(*,'(A)') would do.
 * buflen : length of buffer.
 *
 * Returns
 * -------
 *  1 : a line has been read and stored in buffer.
 *  0 : readline has not been used, either because the package was compiled
 *      without readline support or because stdin or stdout is not a
 *      terminal (e.g. input redirected from a file in a script). The caller
 *      must then read the line with a normal Fortran READ.
 * -1 : end of file (Ctrl-D on an empty line).
 */
int freadline(const char *prompt, char *buffer, int buflen)
{
#ifdef HAVE_READLINE
    static int initialized = 0;
    char *line;
    size_t n;

    if (!isatty(STDIN_FILENO) || !isatty(STDOUT_FILENO))
        return 0;

    if (!initialized) {
        /* allows conditional settings in ~/.inputrc: $if xnirspec ... */
        rl_readline_name = "xnirspec";
        initialized = 1;
    }

    /* default completion function in readline: file names */
    line = readline(prompt);
    if (line == NULL)
        return -1;

    n = strlen(line);
    if (n > (size_t) buflen)
        n = (size_t) buflen;
    memcpy(buffer, line, n);
    memset(buffer + n, ' ', (size_t) buflen - n);

#ifdef HAVE_READLINE_HISTORY_H
    if (line[0] != '\0')
        add_history(line);
#endif
    free(line);
    fflush(stdout);
    return 1;
#else
    (void) prompt;
    (void) buffer;
    (void) buflen;
    return 0;
#endif
}
