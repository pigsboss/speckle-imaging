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
  <meta http-equiv="Content-Type" content="text/html; charset=utf-8">
  <style>
    #terminal {
      font-family: Monospace,"Courier New";
    }
    table {
      border-width: 1px;
      border-spacing: 1px;
      border-style: outset;
      border-collapse: separate;
      font-family: Monospace,"Courier New";
    }
    table td {
      border-width: 1px;
      padding: 1px;
      border-style: inset;
      font-family: Monospace,"Courier New";
    }
  </style>
  <title>Speckle Imaging WebUI - IMAGESC</title>
</head>
  <body>
  <h1>Speckle Imaging WebUI - IMAGESC</h1>
  <form action="imagesc.php" method="post"><table>
    <tr><td>FITS:<td><select name="fits" id="fits">
      <?php
      $lastline=exec('ls',$output);
      foreach($output as $file) {
        if(strtoupper(substr($file,-4))==='FITS') {
          echo '<option label="'.$file.'" value="'.$_SESSION['cwd'].'/'.$file.'">'.$file.'</option>';
        }
      }
      ?>
    </select>
    <tr><td>frame:<td><input type="text" name="frame" id="frame" size="10"></input>(for multi-frame FITS data)
    <tr><td>title:<td><input type="text" name="title" id="title" size="20"></input>(use "\" to escape special characters)
    <tr><td>scale:<td><input type="text" name="scale" id="scale" size="10"></input>(arcsec per pixel)
    <tr><td>x_min:<td><input type="text" name="xmin" id="xmin" size="10"></input>(pixel)
    <tr><td>x_max:<td><input type="text" name="xmax" id="xmax" size="10"></input>(pixel)
    <tr><td>y_min:<td><input type="text" name="ymin" id="ymin" size="10"></input>(pixel)
    <tr><td>y_max:<td><input type="text" name="ymax" id="ymax" size="10"></input>(pixel)
    <tr><td>c_min:<td><input type="text" name="cmin" id="cmin" size="10"></input>
    <tr><td>c_max:<td><input type="text" name="cmax" id="cmax" size="10"></input>
    <tr><td>x axis label:<td><input type="text" name="xlabel" id="xlabel" size="10"></input>
    <tr><td>y axis label:<td><input type="text" name="ylabel" id="ylabel" size="10"></input>
    <tr><td>colorbar label:<td><input type="text" name="clabel" id="clabel" size="10"></input>
    <tr><td>color scheme:<td><select name="color" id="color">
      <option label="rainbow" value="rainbow">Rainbow</option>
      <option label="heat" value="heat">Heat</option>
      <option label="gray" value="gray">Gray</option>
    </select>
    <tr><td>mapping:<td><select name="map" id="map">
      <option label="linear" value="linear">Linear</option>
      <option label="sqrt" value="sqrt">Square root</option>
      <option label="log" value="log">Logarithm</option>
    </select>
    <tr><td colspan="2"><input type="submit" value="Draw"></input><input type="reset" value="Reset"></input>
  </table></form>
  <div id="terminal">
  <?php
  if(count($_POST) > 0) {
    echo '<hr><h2>Environment</h2>';
    $fits=$_POST['fits'];
    echo 'Data: '.$fits.'<br>';
    if(!chdir(INDEXDIR.'/'.WORKDIR)) {
      die('directory '.INDEXDIR.'/'.WORKDIR.' does not exist.');
    }
    echo 'Current working directory: '.getcwd().'<br>';
    $device=basename($fits).'.png';
    echo 'Device: '.$device.'<br>';
    unset($output);
    $cmd=BINDIR.'/imagesc '.$fits;
    if(strlen(trim($_POST['xmin'])) > 0) {
      $cmd=$cmd.' -xmin='.$_POST['xmin'];
    }
    if(strlen(trim($_POST['xmax'])) > 0) {
      $cmd=$cmd.' -xmax='.$_POST['xmax'];
    }
    if(strlen(trim($_POST['ymin'])) > 0) {
      $cmd=$cmd.' -ymin='.$_POST['ymin'];
    }
    if(strlen(trim($_POST['ymax'])) > 0) {
      $cmd=$cmd.' -ymax='.$_POST['ymax'];
    }
    if(strlen(trim($_POST['cmin'])) > 0) {
      $cmd=$cmd.' -min='.$_POST['cmin'];
    }
    if(strlen(trim($_POST['cmax'])) > 0) {
      $cmd=$cmd.' -max='.$_POST['max'];
    }
    if(strlen(trim($_POST['scale'])) > 0) {
      $cmd=$cmd.' -scale='.$_POST['scale'];
    }
    if(strlen(trim($_POST['xlabel'])) > 0) {
      $cmd=$cmd.' -xlabel='.$_POST['xlabel'];
    }
    if(strlen(trim($_POST['ylabel'])) > 0) {
      $cmd=$cmd.' -ylabel='.$_POST['ylabel'];
    }
    if(strlen(trim($_POST['clabel'])) > 0) {
      $cmd=$cmd.' -clabel='.$_POST['clabel'];
    }
    if(strlen(trim($_POST['title'])) > 0) {
      $cmd=$cmd.' -title='.$_POST['title'];
    }
    if(strlen(trim($_POST['map'])) > 0) {
      $cmd=$cmd.' -map='.$_POST['map'];
    }
    if(strlen(trim($_POST['frame'])) > 0) {
      $cmd=$cmd.' -frame='.$_POST['frame'];
    }
    if(strlen(trim($_POST['color'])) > 0) {
      $cmd=$cmd.' -color='.$_POST['color'];
    }
    $cmd=$cmd.' -device='.$device.'/png';
    $lastline=exec($cmd,$output);
    echo '<hr><h2>Output</h2>';
    foreach($output as $line) {
      echo $line.'<br>';
    }
    echo '<img src="./'.WORKDIR.'/'.$device.'" />';
  }
  ?>
  </div>
  </body>
</html>
