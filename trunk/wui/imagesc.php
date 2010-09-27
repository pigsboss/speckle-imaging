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
  <script type="text/javascript" src="jquery.js"></script>
  <meta http-equiv="Content-Type" content="text/html; charset=utf-8">
  <style>
    body {
    }
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
  <div id="form-load"><form action="imagesc.php?action=load" method="post"><table>
    <tr><td>FITS:<td><select name="fits" id="load-fits">
      <?php
      $lastline=exec('ls',$output);
      foreach($output as $file) {
        if(strtoupper(substr($file,-4))==='FITS') {
          echo '<option id="'.md5($file).'" label="'.$file.'" value="'.$file.'">'.$file.'</option>';
        }
      }
      ?>
    </select><input type="submit" value="Load" id="submit-load"></input>
  </table></form></div>
  <script type="text/javascript">
    $("#load-fits").change(function() {
      $("#form-draw").hide();
    });
  </script>
  <div id="form-draw"><form action="imagesc.php?action=draw" method="post"><table>
    <input type="hidden" name="fits" id="draw-fits">
    <input type="hidden" name="nx" id="nx">
    <input type="hidden" name="ny" id="ny">
    <input type="hidden" name="nz" id="nz">
    <input type="hidden" name="min" id="min">
    <input type="hidden" name="max" id="max">
    <tr><td>frame:<td><input type="text" name="frame" id="frame" size="10"></input> (for multi-frame FITS data)
    <tr><td>title:<td><input type="text" name="title" id="title" size="10"></input> (use "\" to escape special characters)
    <tr><td>scale:<td><input type="text" name="scale" id="scale" size="10"></input> arcsec per pixel
    <tr><td>x_min:<td><input type="text" name="xmin" id="xmin" size="10"></input> pixel
    <tr><td>x_max:<td><input type="text" name="xmax" id="xmax" size="10"></input> pixel
    <tr><td>y_min:<td><input type="text" name="ymin" id="ymin" size="10"></input> pixel
    <tr><td>y_max:<td><input type="text" name="ymax" id="ymax" size="10"></input> pixel
    <tr><td>c_min:<td><input type="text" name="cmin" id="cmin" size="10"></input>
    <tr><td>c_max:<td><input type="text" name="cmax" id="cmax" size="10"></input>
    <tr><td>x axis label:<td><input type="text" name="xlabel" id="xlabel" size="10"></input> (default: x/arcsec)
    <tr><td>y axis label:<td><input type="text" name="ylabel" id="ylabel" size="10"></input> (default: y/arcsec)
    <tr><td>colorbar label:<td><input type="text" name="clabel" id="clabel" size="10"></input> (default: Intensity)
    <tr><td>color scheme:<td><select name="color" id="color">
      <option id="colormap-rainbow" label="rainbow" value="rainbow">Rainbow</option>
      <option id="colormap-heat" label="heat" value="heat">Heat</option>
      <option id="colormap-gray" label="gray" value="gray">Gray</option>
    </select>
    <tr><td>mapping:<td><select name="map" id="map">
      <option id="map-linear" label="linear" value="linear">Linear</option>
      <option id="map-sqrt" label="sqrt" value="sqrt">Square root</option>
      <option id="map-log" label="log" value="log">Logarithm</option>
    </select>
    <tr><td colspan="2"><input type="submit" value="Draw"></input><input type="reset" value="Reset"></input>
  </table></form></div>
  <script type="text/javascript">
    $("#form-draw").hide();
  </script>
  <div id="terminal">
  <?php
  if(count($_POST) > 0) {
    if($_REQUEST['action']=='load') {
      echo '<script type="text/javascript">$("#form-draw").show();</script>';
      $fits=$_SESSION['cwd'].'/'.$_POST['fits'];
      echo '<script type="text/javascript">$("#'.md5($_POST['fits']).'").attr("selected","selected");</script>';
      $cmd=BINDIR.'/fitsinfo '.$fits;
      unset($output);
      $lastline=exec($cmd,$output);
      foreach($output as $line) {
        $info=explode(':',$line);
        if(strpos($info[0],'x axis length')) {
          $nx=trim($info[1]);
        }
        if(strpos($info[0],'y axis length')) {
          $ny=trim($info[1]);
        }
        if(strpos($info[0],'z axis length')) {
          $nz=trim($info[1]);
        }
        if(strpos($info[0],'min')) {
          $min=trim($info[1]);
        }
        if(strpos($info[0],'max')) {
          $max=trim($info[1]);
        }
      }
      echo '<script type="text/javascript">$("#frame").after(" (1 - '.$nz.')");</script>';
      echo '<script type="text/javascript">$("#xmin").after(" (1 - '.$nx.')");</script>';
      echo '<script type="text/javascript">$("#xmax").after(" (1 - '.$nx.')");</script>';
      echo '<script type="text/javascript">$("#ymin").after(" (1 - '.$ny.')");</script>';
      echo '<script type="text/javascript">$("#ymax").after(" (1 - '.$ny.')");</script>';
      echo '<script type="text/javascript">$("#cmin").after(" (>= '.$min.')");</script>';
      echo '<script type="text/javascript">$("#cmax").after(" (<= '.$max.')");</script>';
      echo '<script type="text/javascript">$("#submit-load").attr("value","Reload");</script>';
      echo '<script type="text/javascript">$("#draw-fits").attr("value","'.$_POST['fits'].'");</script>';
      echo '<script type="text/javascript">$("#nx").attr("value","'.$nx.'");</script>';
      echo '<script type="text/javascript">$("#ny").attr("value","'.$ny.'");</script>';
      echo '<script type="text/javascript">$("#nz").attr("value","'.$nz.'");</script>';
      echo '<script type="text/javascript">$("#min").attr("value","'.$min.'");</script>';
      echo '<script type="text/javascript">$("#max").attr("value","'.$max.'");</script>';
    }
    if($_REQUEST['action']=='draw') {
      echo '<script type="text/javascript">$("#draw-fits").attr("value","'.$_POST['fits'].'");</script>';
      echo '<script type="text/javascript">$("#nx").attr("value","'.$_POST['nx'].'");</script>';
      echo '<script type="text/javascript">$("#ny").attr("value","'.$_POST['ny'].'");</script>';
      echo '<script type="text/javascript">$("#nz").attr("value","'.$_POST['nz'].'");</script>';
      echo '<script type="text/javascript">$("#min").attr("value","'.$_POST['min'].'");</script>';
      echo '<script type="text/javascript">$("#max").attr("value","'.$_POST['max'].'");</script>';
      echo '<script type="text/javascript">$("#frame").after(" (1 - '.$_POST['nz'].')");</script>';
      echo '<script type="text/javascript">$("#xmin").after(" (1 - '.$_POST['nx'].')");</script>';
      echo '<script type="text/javascript">$("#xmax").after(" (1 - '.$_POST['nx'].')");</script>';
      echo '<script type="text/javascript">$("#ymin").after(" (1 - '.$_POST['ny'].')");</script>';
      echo '<script type="text/javascript">$("#ymax").after(" (1 - '.$_POST['ny'].')");</script>';
      echo '<script type="text/javascript">$("#cmin").after(" (>= '.$_POST['min'].')");</script>';
      echo '<script type="text/javascript">$("#cmax").after(" (<= '.$_POST['max'].')");</script>';
      echo '<script type="text/javascript">$("#form-draw").show();</script>';
      echo '<script type="text/javascript">$("#submit-load").attr("value","Reload");</script>';
      $fits=$_SESSION['cwd'].'/'.$_POST['fits'];
      echo '<script type="text/javascript">$("#'.md5($_POST['fits']).'").attr("selected","selected");</script>';
      echo '<hr><h2>Environment</h2>';
      echo 'Data: '.$fits.'<br>';
      if(!chdir(INDEXDIR.'/'.WORKDIR)) {
        die('directory '.INDEXDIR.'/'.WORKDIR.' does not exist.');
      }
      echo 'Current working directory: '.getcwd().'<br>';
      $cmd=BINDIR.'/imagesc '.$fits;
      if(strlen(trim($_POST['xmin'])) > 0) {
        $cmd=$cmd.' -xmin='.$_POST['xmin'];
        echo '<script type="text/javascript">$("#xmin").attr("value","'.$_POST['xmin'].'");</script>';
      }
      if(strlen(trim($_POST['xmax'])) > 0) {
        $cmd=$cmd.' -xmax='.$_POST['xmax'];
        echo '<script type="text/javascript">$("#xmax").attr("value","'.$_POST['xmax'].'");</script>';
      }
      if(strlen(trim($_POST['ymin'])) > 0) {
        $cmd=$cmd.' -ymin='.$_POST['ymin'];
        echo '<script type="text/javascript">$("#ymin").attr("value","'.$_POST['ymin'].'");</script>';
      }
      if(strlen(trim($_POST['ymax'])) > 0) {
        $cmd=$cmd.' -ymax='.$_POST['ymax'];
        echo '<script type="text/javascript">$("#ymax").attr("value","'.$_POST['ymax'].'");</script>';
      }
      if(strlen(trim($_POST['cmin'])) > 0) {
        $cmd=$cmd.' -min='.$_POST['cmin'];
        echo '<script type="text/javascript">$("#cmin").attr("value","'.$_POST['cmin'].'");</script>';
      }
      if(strlen(trim($_POST['cmax'])) > 0) {
        $cmd=$cmd.' -max='.$_POST['max'];
        echo '<script type="text/javascript">$("#cmax").attr("value","'.$_POST['cmax'].'");</script>';
      }
      if(strlen(trim($_POST['scale'])) > 0) {
        $cmd=$cmd.' -scale='.$_POST['scale'];
        echo '<script type="text/javascript">$("#scale").attr("value","'.$_POST['scale'].'");</script>';
      }
      if(strlen(trim($_POST['xlabel'])) > 0) {
        $cmd=$cmd.' -xlabel='.$_POST['xlabel'];
        echo '<script type="text/javascript">$("#xlabel").attr("value","'.$_POST['xlabel'].'");</script>';
      }
      if(strlen(trim($_POST['ylabel'])) > 0) {
        $cmd=$cmd.' -ylabel='.$_POST['ylabel'];
        echo '<script type="text/javascript">$("#ylabel").attr("value","'.$_POST['ylabel'].'");</script>';
      }
      if(strlen(trim($_POST['clabel'])) > 0) {
        $cmd=$cmd.' -clabel='.$_POST['clabel'];
        echo '<script type="text/javascript">$("#clabel").attr("value","'.$_POST['clabel'].'");</script>';
      }
      if(strlen(trim($_POST['title'])) > 0) {
        $cmd=$cmd.' -title='.$_POST['title'];
        echo '<script type="text/javascript">$("#title").attr("value","'.$_POST['title'].'");</script>';
      }
      if(strlen(trim($_POST['map'])) > 0) {
        $cmd=$cmd.' -map='.$_POST['map'];
        echo '<script type="text/javascript">$("#map-'.$_POST['map'].'").attr("selected","selected");</script>';
      }
      if(strlen(trim($_POST['frame'])) > 0) {
        $cmd=$cmd.' -frame='.$_POST['frame'];
        echo '<script type="text/javascript">$("#frame").attr("value","'.$_POST['frame'].'");</script>';
      }
      if(strlen(trim($_POST['color'])) > 0) {
        $cmd=$cmd.' -color='.$_POST['color'];
        echo '<script type="text/javascript">$("#colormap-'.$_POST['color'].'").attr("selected","selected");</script>';
      }
      $device=basename($fits).'.';
      unset($output);
      $lastline=exec($cmd.' -device='.$device.'png/png',$output);
      $lastline=exec($cmd.' -device='.$device.'ps/vcps',$output);
      echo '<hr><h2>Output</h2>';
      echo '<div><img src="./'.WORKDIR.'/'.$device.'png" /></div>';
      echo '<div><a target="_blank" href="./'.WORKDIR.'/'.$device.'ps">PostScript format</a><br>';
      $lastline=exec('epstopdf '.$device.'ps',$output);
      echo '<a target="_blank" href="./'.WORKDIR.'/'.$device.'pdf">PDF format</a><br></div>';
    }
  }
  ?>
  </div>
  </body>
</html>
