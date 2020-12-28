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
            }
        }

        stage ("deploy")
        {
            steps
            {
                echo "Deploying . . ."
            }
        }
    }
}