<!DOCTYPE html PUBLIC "-//W3C//DTD XHTML 1.0 Strict//EN" "http://www.w3.org/TR/xhtml1/DTD/xhtml1-strict.dtd">
<html lang="en" xmlns="http://www.w3.org/1999/xhtml" xml:lang="en">
<head>
  <meta http-equiv="Content-Type" content="text/html; charset=utf-8">
  <title>Speckle Imaging WebUI</title>
</head>
  <body>
  <h1>Speckle Imaging WebUI</h1>
  <div><form action="index.php" method="post">
    <input type="text" size="40" name="command" id="command"></input><input type="submit" value="Execute"></input>
  </form></div>
  <div>
  <?php
    if(count($_POST) > 0) {
      $command=$_POST['command'];
      echo '<p>Command: '.$command.'<hr>';
      echo '<p>Output:<br>';
      $lastline=exec($command,$output);
      foreach($output as $line) {
        echo $line.'<br>';
      }
    }
  ?>
  </div>
  </body>
</html>
