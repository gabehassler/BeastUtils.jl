module RunBeast #module for running BEAST

export run_beast,
       check_beast,
       find_beast

const BEAST_JAR = "beast.jar"
const BEAST_HOME = "BEAST_HOME"

function find_beast(;beast_home::String="")
    if isempty(beast_home)
        beast_home = haskey(ENV, BEAST_HOME) ? ENV[BEAST_HOME] : pwd()
    end

    if basename(beast_home) == BEAST_JAR
        return abspath(beast_home)
    end

    path = joinpath(beast_home, BEAST_JAR)
    if isfile(path)
        beast_path = path
    else
        path = joinpath(beast_home, "build", "dist", BEAST_JAR)
        if isfile(path)
            beast_path = path
        else
            error("Could not find beast.jar file. Consider setting the " *
                  "$BEAST_HOME environment variable to the location of your " *
                  "BEAST installation or git repo. You may also specify the " *
                  "optional argument 'beast_jar=path/to/your/beast.jar'")
        end
    end

    return abspath(beast_path)
end

function check_beast(;beast_jar::String = find_beast())
    println("Checking Java installation...")
    run(`java -version`)

    println("\n\nChecking BEAST installation...")
    run(`java -jar $beast_jar -version`)

    println("\n\nJava and BEAST checks succeeded.")
end

function run_beast(xml_path::String;
                    seed::Int = -1, # don't set seed by default
                    overwrite::Bool = false, # don't overwrite log files by default
                    directory::String = pwd(), #don't change working directory by default
                    beast_jar::String = find_beast()
                    )

    old_directory = pwd()

    cd(directory)
    cmds = ["java", "-jar", beast_jar]

    if seed != -1
        push!(cmds, "-seed")
        push!(cmds, string(seed))
    end

    if overwrite
        push!(cmds, "-overwrite")
    end

    push!(cmds, "-fail_threads")

    push!(cmds, xml_path)

    out = run(Cmd(cmds))
    if out.exitcode != 0
        error("BEAST run threw an error. See console for error message.")
    end

    cd(old_directory)
end


end
