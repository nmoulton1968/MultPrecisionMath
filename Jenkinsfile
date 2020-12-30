pipeline
{
    agent any 
    
    stages
    {
        stage ("build")
        {
            steps
            {
                echo "Building . . ."
                sh 'make pi'
                sh 'make e'
            }
        }

        stage ("deploy")
        {
            steps
            {
                echo "Deploying . . ."
                fileOperations([fileCopyOperation(
                    excludes: '',
                    flattenFiles: false,
                    includes: '*.exe',
                    targetLocation: "/home/nmoulton/Documents/JenkinsBuilds"
                )])

                echo "Contents of output folder:"
                sh 'ls /home/nmoulton/Documents/JenkinsBuilds/'
            }
        }
    }
}