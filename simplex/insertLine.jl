function skiplines(io::IO, n)
    i = 1
    while i <= n
       eof(io) && error("File contains less than $n lines")
       i += read(io, Char) === '\n'
    end
end

function insertLine(file::String, string::String, lineNr::Number)
    f = open(file, "r+");
    skiplines(f, lineNr);
    skip(f, -1)
    mark(f)
    buf = IOBuffer()
    write(buf, f)
    seekstart(buf)
    reset(f)
    print(f, string);
    write(f, buf)
    close(f)
end