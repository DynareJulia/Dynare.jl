using Pkg
using Pkg.Artifacts

# This is the path to the Artifacts.toml we will manipulate
artifact_toml = joinpath(@__DIR__, "Artifacts.toml")

# Query the `Artifacts.toml` file for the hash bound to the name "dynare-preprocessor"
# (returns `nothing` if no such binding exists)
preprocessor_hash = artifact_hash("dynare-preprocessor", artifact_toml)
@show preprocessor_hash
if preprocessor_hash == nothing || !artifact_exists(preprocessor_hash)
    # create_artifact() returns the content-hash of the artifact directory once we're finished creating it
    preprocessor_hash = create_artifact() do artifact_dir
        @show artifact_dir
        # We create the artifact by simply downloading a few files into the new artifact directory
        preprocessor_url_base = "https://git.dynare.org/julia-packages//-/tree/master/"
        download("$(preprocessor_url_base)/linux/64/preprocess.tar.gz", joinpath(artifact_dir, "preprocessor.tar.gz"))
    end

    # Now bind that hash within our `Artifacts.toml`.  `force = true` means that if it already exists,
    # just overwrite with the new content-hash.  Unless the source files change, we do not expect
    # the content hash to change, so this should not cause unnecessary version control churn.
    bind_artifact!(artifact_toml, "dynare-preprocessor", preprocessor_hash)
end

# Get the path of the iris dataset, either newly created or previously generated.
# this should be something like `~/.julia/artifacts/dbd04e28be047a54fbe9bf67e934be5b5e0d357a`
dynare_preprocessor_path = artifact_path(preprocessor_hash)
@show dynare_preprocessor_path
