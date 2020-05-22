module RunBeast #module for running BEAST

const BEAST_JAR = "beast.jar"

function find_beast(beast_home::String)
    path = joinpath(beast_home, BEAST_JAR)
    if isfile(path)
        beast_path = path
    else
        path = joinpath(beast_home, "build", "dist", BEAST_JAR)
        if isfile(path)
            beast_path = path
        else
            error("Could not find BEAST .jar file.") # TODO: better error message
        end
    end

    return beast_path
end

function check_beast(;beast_jar::String = find_beast(ENV["BEAST_HOME"]))
    println("Checking Java installation...")
    run(`java -version`)

    println("Checking BEAST installation...")
    run(`java -jar $beast_jar -version`)

    println("Java and BEAST checks suceeded.")
end

function run_beast(xml_path::String;
                    seed::Int = -1, # don't set seed by default
                    overwrite::Bool = false, # don't overwrite log files by default
                    directory::String = pwd(), #don't change working directory by default
                    beast_jar::String = find_beast(ENV["BEAST_HOME"])
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

    push!(cmds, xml_path)

    run(Cmd(cmds))

    cd(old_directory)

end



end