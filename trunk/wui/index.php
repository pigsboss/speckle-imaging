<?php
session_start();
require_once('config.php');
if(array_key_exists('cwd',$_SESSION)) {
  if(!chdir($_SESSION['cwd'])) {
    die('directory '.$_SESSION['cwd'].' does not exist.');
  }
}else {
  if(!chdir(INDEXDIR.'/'.WORKDIR)) {
    die('directory '.INDEXDIR.'/'.WORKDIR.' does not exist.');
  }
}
?>
<!DOCTYPE html PUBLIC "-//W3C//DTD XHTML 1.0 Strict//EN" "http://www.w3.org/TR/xhtml1/DTD/xhtml1-strict.dtd">
<html lang="en" xmlns="http://www.w3.org/1999/xhtml" xml:lang="en">
<head>
  <script type="text/javascript" src="_lib/jquery.js"></script>
  <script type="text/javascript" src="jquery.jstree.js"></script>
  <meta http-equiv="Content-Type" content="text/html; charset=utf-8">
  <style>
    #terminal {
      font-family: "Monospace","Courier New";
    }
  </style>
  <title>Speckle Imaging WebUI</title>
</head>
  <body>
  <h1>Speckle Imaging WebUI</h1>
  <div><form action="index.php" method="post">
    <input type="text" size="40" name="dir" id="dir"></input><input type="submit" value="Change Working Directory"></input><br>
    <input type="text" size="40" name="command" id="command"></input><input type="submit" value="Execute"></input>
  </form></div>
  <div>
    <input type="button" value="IMAGESC" onclick="window.open('imagesc.php','IMAGESC','scrollbars=yes,toolbar=no,resizeable=yes,directories=no,menubar=no,location=no');"></input>  
    </div>
  <div id="terminal">
  <?php
  if(count($_POST) > 0) {
    $command=$_POST['command'];
    $dir=$_POST['dir'];
    if(strlen(trim($dir)) > 0) {
      if(!chdir($dir)) {
        die('directory '.$dir.' does not exist.');
      }
      $_SESSION['cwd']=getcwd();
    }else if(strlen(trim($command)) > 0) {
      echo '<hr><h2>Command</h2>'.$command.'<hr>';
      echo '<h2>Output</h2>';
      unset($output);
      $lastline=exec($command,$output);
      foreach($output as $line) {
        echo $line.'<br>';
      }
    }
  }
  ?>
  </div>
  </body>
</html>
